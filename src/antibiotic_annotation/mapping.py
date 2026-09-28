"""Chemical identity -> ChEBI id(s) via UniChem.

Order of preference
1. ``direct_ccd_crossref``: a ChEBI cross-link on the compound record, accepted only when that
   ChEBI entry's standard InChIKey equals the deposited one.
2. ``unichem_inchikey``: exact standard InChIKey match. Several ChEBI ids can share one standard
   InChIKey (tautomers, zwitterions, polymer classes); they are treated as the same *standardised
   chemical identity for mapping purposes*: a deterministic primary is chosen and ontology
   evidence is unioned.
3. ``unichem_connectivity_protonation``: UniChem connectivity search, accepted only when the
   difference is limited to protonation / charge / H-atom count.
4. ``unichem_connectivity_stereo_undefined``: the deposited descriptor has no stereo layer
   (typical for BIRD entries) and exactly one plausible single-species ChEBI identity matches.
5. ``manual_synonym``: curated table, used only when no structural identifier exists.

Stereo or connectivity conflicts are never resolved automatically: status
``unresolved_identity_conflict`` with every candidate and its mismatch flags retained.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from .cache import HttpError, JsonFileCache, NetworkUnavailable
from .models import (
    ENTITY_CHEBI_INPUT,
    MAP_RESOLVED,
    MAP_UNRESOLVED_CONFLICT,
    MAP_UNRESOLVED_NETWORK,
    MAP_UNRESOLVED_NO_CHEBI,
    MAP_UNRESOLVED_NO_STRUCTURE,
    METHOD_CHEBI_INPUT,
    METHOD_DIRECT_XREF,
    METHOD_MANUAL_SYNONYM,
    METHOD_UNICHEM_CONNECTIVITY,
    METHOD_UNICHEM_INCHIKEY,
    METHOD_UNICHEM_STEREO_UNDEFINED,
    ChemicalIdentity,
    MappingCandidate,
    MappingResult,
    normalise_chebi_id,
)

UNICHEM_COMPOUNDS_URL = "https://www.ebi.ac.uk/unichem/api/v1/compounds"
UNICHEM_CONNECTIVITY_URL = "https://www.ebi.ac.uk/unichem/api/v1/connectivity"

# UniChem comparison flags that may be False for a "same compound, different protonation state /
# salt form" match. Anything else (connectivity, stereo, isotope, formula) is an identity conflict.
ACCEPTABLE_MISMATCHES = frozenset({"protonation", "charge", "HAtoms", "isotopicExchangeableH"})
STEREO_MISMATCHES = frozenset({"stereoSp3", "stereoSp3Inverted", "stereoType", "stereoDbond"})
NO_STEREO_LAYER = "UHFFFAOYSA"  # second InChIKey block when no stereo is defined

CONFIDENCE = {
    METHOD_DIRECT_XREF: "high",
    METHOD_UNICHEM_INCHIKEY: "high",
    METHOD_CHEBI_INPUT: "high",
    METHOD_UNICHEM_CONNECTIVITY: "medium",
    METHOD_UNICHEM_STEREO_UNDEFINED: "medium",
    METHOD_MANUAL_SYNONYM: "low",
}
RANK_HOW = "most specific entity (not a macromolecule/mixture class, fewest descendants), highest star rating, then lowest numeric id"


@dataclass
class ChebiLookup:
    """What the mapper needs to know about a ChEBI id (supplied by the ontology layer)."""

    stars: int | None = None
    inchikey: str | None = None
    is_class_like: bool = False          # macromolecule / polymer / mixture class, not a single species
    num_descendants: int | None = None
    forms: frozenset[str] = frozenset()  # ids linked by conjugate acid/base or tautomer relations


class UniChemClient:
    def __init__(self, transport, cache: JsonFileCache):
        self.transport = transport
        self.cache = cache

    def exact(self, inchikey: str) -> dict[str, Any]:
        return self.cache.get_or_fetch("unichem/inchikey", inchikey, lambda: self.transport.post_json(UNICHEM_COMPOUNDS_URL, {"type": "inchikey", "compound": inchikey}))

    def connectivity(self, inchikey: str) -> dict[str, Any]:
        payload = {"type": "inchikey", "compound": inchikey, "searchComponents": True, "searchSources": True}
        return self.cache.get_or_fetch("unichem/connectivity", inchikey, lambda: self.transport.post_json(UNICHEM_CONNECTIVITY_URL, payload))


def _num(chebi_id: str) -> int:
    return int(chebi_id.split(":")[1])


def stereo_undefined(inchikey: str) -> bool:
    parts = inchikey.split("-")
    return len(parts) > 1 and parts[1].startswith(NO_STEREO_LAYER)


def _chebi_ids_from_exact(raw: dict[str, Any]) -> list[str]:
    ids: list[str] = []
    for comp in raw.get("compounds") or []:
        for s in comp.get("sources") or []:
            cid = normalise_chebi_id(s.get("compoundId")) if s.get("shortName") == "chebi" else None
            if cid and cid not in ids:
                ids.append(cid)
    return ids


def _chebi_candidates_from_connectivity(raw: dict[str, Any]) -> list[MappingCandidate]:
    """One candidate per ChEBI id, keeping the closest match when UniChem lists several."""
    best: dict[str, MappingCandidate] = {}
    for s in raw.get("sources") or []:
        cid = normalise_chebi_id(s.get("compoundId")) if s.get("shortName") == "chebi" else None
        if not cid:
            continue
        comparison = {k: bool(v) for k, v in (s.get("comparison") or {}).items()}
        mismatches = sorted(k for k, v in comparison.items() if not v)
        cand = MappingCandidate(chebi_id=cid, inchikey=s.get("inchikey"), comparison=comparison, mismatches=mismatches, accepted=set(mismatches) <= ACCEPTABLE_MISMATCHES)
        if cid not in best or len(mismatches) < len(best[cid].mismatches):
            best[cid] = cand
    return sorted(best.values(), key=lambda c: (len(c.mismatches), _num(c.chebi_id)))


def load_synonym_table(path: str | Path | None) -> dict[str, dict[str, str]]:
    if path is None or not Path(path).exists():
        return {}
    return json.loads(Path(path).read_text(encoding="utf-8"))


class Mapper:
    """Maps a :class:`ChemicalIdentity` to ChEBI.

    ``lookup`` (ChEBI id -> :class:`ChebiLookup`) is supplied by the ontology layer and is used to
    verify cross-links, rank equivalent ids and group protonation variants of one entity.
    """

    def __init__(self, unichem: UniChemClient, lookup: Callable[[str], ChebiLookup] | None = None, synonym_table: dict[str, dict[str, str]] | None = None):
        self.unichem = unichem
        self.lookup = lookup
        self.synonyms = synonym_table or {}

    # -- helpers ---------------------------------------------------------------
    def _info(self, cid: str) -> ChebiLookup:
        return self.lookup(cid) if self.lookup else ChebiLookup()

    def _rank(self, cid: str) -> tuple:
        i = self._info(cid)
        return (int(i.is_class_like), i.num_descendants if i.num_descendants is not None else 10**6, -(i.stars or 0), _num(cid))

    def _resolved(self, method: str, source_id: str, ids: list[str], candidates: list[MappingCandidate], note: str, union_note: str, primary: str | None = None) -> MappingResult:
        """Build a resolved result over ids that share one standardised identity."""
        how = RANK_HOW if self.lookup else "lowest numeric id"
        if primary is None:
            primary = min(ids, key=self._rank)
        else:
            how = "fewest InChIKey-layer differences, then " + how
        res = MappingResult(status=MAP_RESOLVED, method=method, source_id=source_id, primary_chebi_id=primary, equivalent_chebi_ids=ids, candidates=candidates, confidence=CONFIDENCE[method])
        res.notes.append(f"{note}; primary chosen by {how}")
        if len(ids) > 1:
            res.evidence_unioned = True
            res.notes.append(union_note)
        return res

    def _conflict(self, ident: ChemicalIdentity, cands: list[MappingCandidate], extra_notes: list[str]) -> MappingResult:
        res = MappingResult(status=MAP_UNRESOLVED_CONFLICT, source_id=ident.inchikey, candidates=cands, notes=list(extra_notes))
        res.notes.append("ChEBI entries share the InChIKey connectivity layer but differ in " + "; ".join(f"{c.chebi_id}: {','.join(c.mismatches)}" for c in cands) + " -- not mapped automatically")
        return self._synonym_fallback(ident, res, allow=False)

    def _group_forms(self, ids: list[str]) -> list[set[str]]:
        """Union-find over conjugate acid/base and tautomer links: one group per chemical identity."""
        groups: list[set[str]] = []
        for cid in ids:
            linked = self._info(cid).forms | {cid}
            merged, rest = {cid}, []
            for g in groups:
                if g & linked:
                    merged |= g
                else:
                    rest.append(g)
            groups = rest + [merged]
        return groups

    def _synonym_fallback(self, ident: ChemicalIdentity, res: MappingResult, allow: bool = True) -> MappingResult:
        """Curated synonym mapping: only when no structural identifier exists (``allow``)."""
        entry = self.synonyms.get(ident.input_id.upper())
        cid = normalise_chebi_id(entry.get("chebi_id")) if entry else None
        if cid and allow:
            return MappingResult(status=MAP_RESOLVED, method=METHOD_MANUAL_SYNONYM, source_id=ident.input_id, primary_chebi_id=cid, equivalent_chebi_ids=[cid], confidence=CONFIDENCE[METHOD_MANUAL_SYNONYM], notes=[f"manual synonym table: {entry.get('note', '')}"])
        if entry and not allow:
            res.notes.append("a manual synonym entry exists but is ignored because a structural identifier is available")
        return res

    # -- main ------------------------------------------------------------------
    def map(self, ident: ChemicalIdentity) -> MappingResult:
        if ident.entity_kind == ENTITY_CHEBI_INPUT:
            cid = normalise_chebi_id(ident.input_id)
            return MappingResult(status=MAP_RESOLVED, method=METHOD_CHEBI_INPUT, source_id=cid, primary_chebi_id=cid, equivalent_chebi_ids=[cid], confidence=CONFIDENCE[METHOD_CHEBI_INPUT])

        # 1. ChEBI cross-links on the compound record, verified by InChIKey equality (database
        #    cross-links alone are not identity-verified: e.g. water linked to "oxygen atom").
        xref_candidates: list[MappingCandidate] = []
        for cid in filter(None, map(normalise_chebi_id, ident.xrefs.get("ChEBI", []))):
            key = self._info(cid).inchikey
            if ident.inchikey and key == ident.inchikey:
                mismatches, why = [], "verified by standard InChIKey equality"
            else:
                mismatches = ["inchikey_not_equal" if key else "chebi_entry_has_no_structure"]
                why = f"NOT verified ({mismatches[0]})"
            xref_candidates.append(MappingCandidate(chebi_id=cid, inchikey=key, accepted=not mismatches, mismatches=mismatches, note=f"ChEBI cross-reference on the compound record, {why}"))
        verified = [c.chebi_id for c in xref_candidates if c.accepted]
        if verified:
            return self._resolved(METHOD_DIRECT_XREF, ident.input_id, verified, xref_candidates, "compound-record ChEBI cross-reference verified by InChIKey equality", "several verified entries share the standard InChIKey: same standardised chemical identity for mapping purposes; ontology evidence unioned")

        if not ident.inchikey:
            return self._synonym_fallback(ident, MappingResult(status=MAP_UNRESOLVED_NO_STRUCTURE, source_id=ident.input_id, notes=["no InChIKey available for the deposited species"]))

        # 2. exact standard InChIKey
        try:
            ids = _chebi_ids_from_exact(self.unichem.exact(ident.inchikey))
        except (NetworkUnavailable, HttpError) as exc:
            return MappingResult(status=MAP_UNRESOLVED_NETWORK, source_id=ident.inchikey, notes=[f"UniChem exact lookup failed: {exc}"])
        if ids:
            cands = [MappingCandidate(chebi_id=i, inchikey=ident.inchikey, accepted=True, note="exact standard InChIKey match") for i in ids] + xref_candidates
            return self._resolved(METHOD_UNICHEM_INCHIKEY, ident.inchikey, ids, cands, "exact standard InChIKey match",
                                  f"{len(ids)} ChEBI entries share standard InChIKey {ident.inchikey}: treated as the same standardised chemical identity for mapping purposes, ontology evidence unioned over all of them")

        # 3. connectivity search (salt / protonation forms; stereo-undefined descriptors)
        try:
            cands = _chebi_candidates_from_connectivity(self.unichem.connectivity(ident.inchikey))
        except (NetworkUnavailable, HttpError) as exc:
            return MappingResult(status=MAP_UNRESOLVED_NETWORK, source_id=ident.inchikey, notes=[f"UniChem connectivity lookup failed: {exc}"])
        no_stereo = stereo_undefined(ident.inchikey)
        known = {c.chebi_id for c in cands}
        cands += [x for x in xref_candidates if x.chebi_id not in known]
        accepted = [c for c in cands if c.accepted]
        if accepted:
            primary = min(accepted, key=lambda c: (len(c.mismatches), self._rank(c.chebi_id))).chebi_id
            diffs = sorted({m for c in accepted for m in c.mismatches})
            return self._resolved(METHOD_UNICHEM_CONNECTIVITY, ident.inchikey, [c.chebi_id for c in accepted], cands,
                                  f"no exact InChIKey match; accepted connectivity match differing only in {diffs or ['nothing']}", "evidence unioned over accepted protonation/charge variants", primary=primary)
        if not cands:
            return self._synonym_fallback(ident, MappingResult(status=MAP_UNRESOLVED_NO_CHEBI, source_id=ident.inchikey, notes=["no ChEBI entry found by exact InChIKey or connectivity search"]), allow=False)
        if not no_stereo:
            return self._conflict(ident, cands, [])

        # 4. stereo cannot be compared: accept only a single plausible single-species identity,
        #    counting protonation variants of one entity as one group.
        for c in cands:
            c.note = "deposited descriptor has no stereo layer (UHFFFAOYSA): stereo cannot be compared"
        plausible = [c.chebi_id for c in cands if c.mismatches and set(c.mismatches) <= STEREO_MISMATCHES | ACCEPTABLE_MISMATCHES and not self._info(c.chebi_id).is_class_like]
        groups = self._group_forms(plausible)
        if len(groups) == 1:
            ids = sorted(groups[0], key=self._rank)
            for c in cands:
                c.accepted = c.chebi_id in ids
                if c.accepted:
                    c.note = "stereo-undefined query: accepted as the single plausible connectivity candidate"
            return self._resolved(METHOD_UNICHEM_STEREO_UNDEFINED, ident.inchikey, ids, cands, "deposited descriptor has no stereo layer; connectivity matches exactly one plausible non-class ChEBI identity; NOT equivalent to an exact InChIKey mapping", "evidence unioned over protonation variants of that identity")
        notes = ["conflict_reason=stereo_undefined_in_query"]
        if len(groups) > 1:
            notes.append(f"multiple_stereo_candidates={len(groups)}: " + " vs ".join("/".join(sorted(g)) for g in groups))
        return self._conflict(ident, cands, notes)
