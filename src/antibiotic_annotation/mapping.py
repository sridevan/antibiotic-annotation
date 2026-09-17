"""Phase 3, step 2: chemical identity -> ChEBI id(s) via UniChem.

Order of preference
1. ``direct_ccd_crossref``: a ChEBI cross-reference asserted by RCSB on the CCD (rare).
2. ``unichem_inchikey``: exact standard InChIKey match. Several ChEBI ids can share one standard
   InChIKey (tautomers, zwitterions); they are treated as the same *standardised chemical identity
   for mapping purposes*, a deterministic primary is chosen and ontology evidence is unioned.
3. ``unichem_connectivity_protonation``: UniChem connectivity search accepted only when the
   difference is limited to protonation / charge / H-atom count. Stereo or connectivity
   differences are never accepted automatically: status ``unresolved_identity_conflict`` with the
   candidates and their mismatch flags retained.
4. ``manual_synonym``: curated table, used only when no structural identifier exists.
"""
from __future__ import annotations

import json
from pathlib import Path
from dataclasses import dataclass
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
CACHE_NS_EXACT = "unichem/inchikey"
CACHE_NS_CONN = "unichem/connectivity"

# UniChem comparison flags that may be False for an accepted "same compound, different
# protonation state / salt form" match. Anything else (connectivity, stereo, isotope, formula)
# is an identity conflict.
ACCEPTABLE_MISMATCHES = frozenset({"protonation", "charge", "HAtoms", "isotopicExchangeableH"})
# UniChem flags that concern only the stereo layer of the InChIKey
STEREO_MISMATCHES = frozenset({"stereoSp3", "stereoSp3Inverted", "stereoType", "stereoDbond"})

CONFIDENCE = {
    METHOD_DIRECT_XREF: "high",
    METHOD_UNICHEM_INCHIKEY: "high",
    METHOD_CHEBI_INPUT: "high",
    METHOD_UNICHEM_CONNECTIVITY: "medium",
    METHOD_UNICHEM_STEREO_UNDEFINED: "medium",
    METHOD_MANUAL_SYNONYM: "low",
}


@dataclass
class ChebiLookup:
    """What the mapper needs to know about a ChEBI id (supplied by the ontology layer)."""

    stars: int | None = None
    inchikey: str | None = None
    is_class_like: bool = False       # macromolecule / polymer / mixture class, not a single species
    num_descendants: int | None = None
    forms: frozenset[str] = frozenset()  # ids linked by conjugate acid/base or tautomer relations


class UniChemClient:
    def __init__(self, transport, cache: JsonFileCache):
        self.transport = transport
        self.cache = cache

    def exact(self, inchikey: str) -> dict[str, Any]:
        return self.cache.get_or_fetch(
            CACHE_NS_EXACT, inchikey,
            lambda: self.transport.post_json(UNICHEM_COMPOUNDS_URL, {"type": "inchikey", "compound": inchikey}),
        )

    def connectivity(self, inchikey: str) -> dict[str, Any]:
        return self.cache.get_or_fetch(
            CACHE_NS_CONN, inchikey,
            lambda: self.transport.post_json(
                UNICHEM_CONNECTIVITY_URL,
                {"type": "inchikey", "compound": inchikey, "searchComponents": True, "searchSources": True},
            ),
        )


def _chebi_ids_from_exact(raw: dict[str, Any]) -> list[tuple[str, str | None]]:
    out: list[tuple[str, str | None]] = []
    for comp in raw.get("compounds") or []:
        key = comp.get("standardInchiKey") or comp.get("inchikey")
        for s in comp.get("sources") or []:
            if s.get("shortName") == "chebi":
                cid = normalise_chebi_id(s.get("compoundId"))
                if cid and cid not in [o[0] for o in out]:
                    out.append((cid, key))
    return out


def _chebi_candidates_from_connectivity(raw: dict[str, Any]) -> list[MappingCandidate]:
    cands: list[MappingCandidate] = []
    for s in raw.get("sources") or []:
        if s.get("shortName") != "chebi":
            continue
        cid = normalise_chebi_id(s.get("compoundId"))
        if not cid:
            continue
        comparison = {k: bool(v) for k, v in (s.get("comparison") or {}).items()}
        mismatches = sorted(k for k, v in comparison.items() if not v)
        key = None
        for k in ("inchikey", "standardInchiKey", "inchiKey"):
            if s.get(k):
                key = s[k]
        if not key and isinstance(s.get("compounds"), list) and s["compounds"]:
            key = s["compounds"][0].get("inchikey") or s["compounds"][0].get("standardInchiKey")
        cand = MappingCandidate(chebi_id=cid, inchikey=key, comparison=comparison, mismatches=mismatches)
        cand.accepted = bool(mismatches) and set(mismatches) <= ACCEPTABLE_MISMATCHES or (not mismatches)
        cands.append(cand)
    # deduplicate by chebi id, keeping the closest match (fewest mismatches)
    best: dict[str, MappingCandidate] = {}
    for c in cands:
        if c.chebi_id not in best or len(c.mismatches) < len(best[c.chebi_id].mismatches):
            best[c.chebi_id] = c
    return sorted(best.values(), key=lambda c: (len(c.mismatches), _num(c.chebi_id)))


def _num(chebi_id: str) -> int:
    return int(chebi_id.split(":")[1])


def load_synonym_table(path: str | Path | None) -> dict[str, dict[str, str]]:
    if path is None or not Path(path).exists():
        return {}
    with Path(path).open(encoding="utf-8") as fh:
        return json.load(fh)


class Mapper:
    """Maps a :class:`ChemicalIdentity` to ChEBI. ``rank`` optionally supplies (stars) per ChEBI id
    so the deterministic primary among equivalent ids prefers better-curated entries."""

    def __init__(self, unichem: UniChemClient, lookup: Callable[[str], ChebiLookup] | None = None, synonym_table: dict[str, dict[str, str]] | None = None):
        self.unichem = unichem
        self.lookup = lookup  # ChEBI id -> ChebiLookup (verification of cross-references, ranking, form grouping)
        self.synonyms = synonym_table or {}

    def _info(self, cid: str) -> ChebiLookup:
        return self.lookup(cid) if self.lookup else ChebiLookup()

    def _rank(self, cid: str) -> tuple:
        i = self._info(cid)
        return (int(i.is_class_like), i.num_descendants if i.num_descendants is not None else 10**6, -(i.stars or 0), _num(cid))

    @staticmethod
    def _finish(res: MappingResult) -> MappingResult:
        res.confidence = CONFIDENCE.get(res.method, "none") if res.resolved else "none"
        return res

    # -- helpers ---------------------------------------------------------------
    def _choose_primary(self, ids: list[str]) -> tuple[str, str]:
        """Deterministic primary among ids sharing one standardised identity.

        With a ranker: most specific entity first (not a macromolecule/mixture class, fewest
        is_a descendants), then highest star rating, then lowest numeric id. Without: lowest id.
        """
        ordered = sorted(ids, key=self._rank)
        how = "most specific entity (not a macromolecule/mixture class, fewest descendants), highest star rating, then lowest numeric id" if self.lookup else "lowest numeric id"
        return ordered[0], how

    # -- main ------------------------------------------------------------------
    def map(self, ident: ChemicalIdentity) -> MappingResult:
        return self._finish(self._map(ident))

    def _map(self, ident: ChemicalIdentity) -> MappingResult:
        res = MappingResult()
        if ident.entity_kind == ENTITY_CHEBI_INPUT:
            cid = normalise_chebi_id(ident.input_id)
            return MappingResult(status=MAP_RESOLVED, method=METHOD_CHEBI_INPUT, source_id=cid, primary_chebi_id=cid, equivalent_chebi_ids=[cid])

        # 1. cross-references asserted on the CCD by RCSB. These are "assigned by PubChem
        #    resource" and are NOT identity-verified (water links to "oxygen atom", glucose to
        #    glucans), so a cross-reference is accepted only when the ChEBI entry's standard
        #    InChIKey equals the deposited species' InChIKey.
        xref_ids = [x for x in (normalise_chebi_id(x) for x in ident.xrefs.get("ChEBI", [])) if x]
        xref_candidates: list[MappingCandidate] = []
        verified: list[str] = []
        for cid in xref_ids:
            key = self._info(cid).inchikey
            cand = MappingCandidate(chebi_id=cid, inchikey=key, note="RCSB cross-reference (assigned by PubChem resource)")
            if ident.inchikey and key == ident.inchikey:
                cand.accepted = True
                cand.note += ", verified by standard InChIKey equality"
                verified.append(cid)
            else:
                cand.accepted = False
                cand.mismatches = ["inchikey_not_equal" if key else "chebi_entry_has_no_structure"]
                cand.note += ", NOT verified (" + cand.mismatches[0] + ")"
            xref_candidates.append(cand)
        if verified:
            primary, how = self._choose_primary(verified)
            res = MappingResult(status=MAP_RESOLVED, method=METHOD_DIRECT_XREF, source_id=ident.input_id, primary_chebi_id=primary, equivalent_chebi_ids=verified, candidates=xref_candidates)
            res.notes.append(f"RCSB ChEBI cross-reference verified by InChIKey equality; primary chosen by {how}")
            if len(verified) > 1:
                res.evidence_unioned = True
                res.notes.append("several verified entries share the standard InChIKey: same standardised chemical identity for mapping purposes; ontology evidence unioned")
            return res

        if not ident.inchikey:
            return self._synonym_fallback(ident, MappingResult(status=MAP_UNRESOLVED_NO_STRUCTURE, source_id=ident.input_id, notes=["no InChIKey available for the deposited species"]))

        # 2. exact standard InChIKey
        try:
            exact_raw = self.unichem.exact(ident.inchikey)
        except (NetworkUnavailable, HttpError) as exc:
            return MappingResult(status=MAP_UNRESOLVED_NETWORK, source_id=ident.inchikey, notes=[f"UniChem exact lookup failed: {exc}"])
        hits = _chebi_ids_from_exact(exact_raw)
        if hits:
            ids = [h[0] for h in hits]
            primary, how = self._choose_primary(ids)
            res = MappingResult(status=MAP_RESOLVED, method=METHOD_UNICHEM_INCHIKEY, source_id=ident.inchikey, primary_chebi_id=primary, equivalent_chebi_ids=ids)
            res.candidates = [MappingCandidate(chebi_id=i, inchikey=ident.inchikey, accepted=True, note="exact standard InChIKey match") for i in ids] + xref_candidates
            if len(ids) > 1:
                res.evidence_unioned = True
                res.notes.append(
                    f"{len(ids)} ChEBI entries share standard InChIKey {ident.inchikey}; treated as the same standardised chemical identity for mapping purposes, primary chosen by {how}, ontology evidence unioned over all of them"
                )
            return res

        # 3. connectivity search (salt / protonation forms)
        try:
            conn_raw = self.unichem.connectivity(ident.inchikey)
        except (NetworkUnavailable, HttpError) as exc:
            return MappingResult(status=MAP_UNRESOLVED_NETWORK, source_id=ident.inchikey, notes=[f"UniChem connectivity lookup failed: {exc}"])
        cands = _chebi_candidates_from_connectivity(conn_raw)
        query_stereo_undefined = ident.inchikey.split("-")[1].startswith("UHFFFAOYSA") if "-" in ident.inchikey else False
        for c in cands:
            if not c.accepted and query_stereo_undefined:
                c.note = "deposited descriptor has no stereo layer (UHFFFAOYSA): stereo cannot be compared; not mapped automatically"
        cands = cands + [x for x in xref_candidates if x.chebi_id not in {c.chebi_id for c in cands}]
        accepted = [c for c in cands if c.accepted]
        if accepted:
            # closest match first (fewest differences), then star rating, then lowest id
            closest = min(len(c.mismatches) for c in accepted)
            ids = [c.chebi_id for c in accepted]
            primary, how = self._choose_primary([c.chebi_id for c in accepted if len(c.mismatches) == closest])
            how = "fewest InChIKey-layer differences, then " + how
            res = MappingResult(status=MAP_RESOLVED, method=METHOD_UNICHEM_CONNECTIVITY, source_id=ident.inchikey, primary_chebi_id=primary, equivalent_chebi_ids=ids, candidates=cands)
            diffs = sorted({m for c in accepted for m in c.mismatches})
            res.notes.append(f"no exact InChIKey match; accepted connectivity match differing only in {diffs or ['nothing']}; primary chosen by {how}")
            if len(ids) > 1:
                res.evidence_unioned = True
                res.notes.append("evidence unioned over accepted protonation/charge variants")
            return res
        if cands and query_stereo_undefined:
            # The deposited descriptor carries no stereo layer, so stereo cannot be compared.
            # Accept only when exactly one plausible identity remains: connectivity matches, the
            # differences are limited to stereo (+ protonation) flags, the entry is a single
            # species (not a class), and protonation variants of one entity count as one group.
            plausible = [c for c in cands if c.mismatches and set(c.mismatches) <= (STEREO_MISMATCHES | ACCEPTABLE_MISMATCHES) and not self._info(c.chebi_id).is_class_like]
            groups = self._group_forms([c.chebi_id for c in plausible])
            if len(groups) == 1:
                ids = sorted(groups[0], key=self._rank)
                for c in cands:
                    c.accepted = c.chebi_id in ids
                    if c.accepted:
                        c.note = "stereo-undefined query: accepted as the single plausible connectivity candidate"
                primary, how = self._choose_primary(ids)
                res = MappingResult(status=MAP_RESOLVED, method=METHOD_UNICHEM_STEREO_UNDEFINED, source_id=ident.inchikey, primary_chebi_id=primary, equivalent_chebi_ids=ids, candidates=cands)
                res.notes.append(f"deposited descriptor has no stereo layer; connectivity matches exactly one plausible non-class ChEBI identity; NOT equivalent to an exact InChIKey mapping; primary chosen by {how}")
                if len(ids) > 1:
                    res.evidence_unioned = True
                    res.notes.append("evidence unioned over protonation variants of that identity")
                return res
            res = MappingResult(status=MAP_UNRESOLVED_CONFLICT, source_id=ident.inchikey, candidates=cands)
            res.notes.append("conflict_reason=stereo_undefined_in_query")
            if len(groups) > 1:
                res.notes.append(f"multiple_stereo_candidates={len(groups)}: " + " vs ".join("/".join(sorted(g)) for g in groups))
            res.notes.append("ChEBI entries share the InChIKey connectivity layer but differ in " + "; ".join(f"{c.chebi_id}: {','.join(c.mismatches)}" for c in cands) + " -- not mapped automatically")
            return self._synonym_fallback(ident, res, allow=False)
        if cands:
            res = MappingResult(status=MAP_UNRESOLVED_CONFLICT, source_id=ident.inchikey, candidates=cands)
            res.notes.append("ChEBI entries share the InChIKey connectivity layer but differ in " + "; ".join(f"{c.chebi_id}: {','.join(c.mismatches)}" for c in cands) + " -- not mapped automatically")
            return self._synonym_fallback(ident, res, allow=False)
        return self._synonym_fallback(ident, MappingResult(status=MAP_UNRESOLVED_NO_CHEBI, source_id=ident.inchikey, candidates=cands, notes=["no ChEBI entry found by exact InChIKey or connectivity search"]), allow=False)

    def _group_forms(self, ids: list[str]) -> list[set[str]]:
        """Union-find over conjugate acid/base and tautomer links: one group per chemical identity."""
        groups: list[set[str]] = []
        for cid in ids:
            linked = self._info(cid).forms | {cid}
            merged = {cid}
            rest = []
            for g in groups:
                if g & linked or cid in g:
                    merged |= g
                else:
                    rest.append(g)
            groups = rest + [merged]
        return groups

    def _synonym_fallback(self, ident: ChemicalIdentity, res: MappingResult, allow: bool = True) -> MappingResult:
        """Curated synonym mapping: only when no structural identifier exists (``allow``)."""
        entry = self.synonyms.get(ident.input_id.upper())
        if entry and allow:
            cid = normalise_chebi_id(entry.get("chebi_id"))
            if cid:
                return MappingResult(status=MAP_RESOLVED, method=METHOD_MANUAL_SYNONYM, source_id=ident.input_id, primary_chebi_id=cid, equivalent_chebi_ids=[cid], notes=[f"manual synonym table: {entry.get('note', '')}"])
        elif entry and not allow:
            res.notes.append("a manual synonym entry exists but is ignored because a structural identifier is available")
        return res
