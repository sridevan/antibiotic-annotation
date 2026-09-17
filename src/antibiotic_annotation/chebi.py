"""Phase 4: ChEBI term retrieval, local cache and graph closures.

Primary source: OLS4 v2 (graph: direct parents, has-role restrictions, precomputed ancestors,
labels, definitions). Secondary source: the ChEBI backend API (star rating, InChIKey, incoming
relations such as "has part" from mixtures, which give family/component context).

Every term is stored once as ``cache/chebi/CHEBI_<n>.json`` (normalised record + trimmed raw
payloads). The is_a closure is computed locally from cached direct parents so that the
benchmark is reproducible offline; OLS4's precomputed ancestor list is kept for cross-checking.
"""
from __future__ import annotations

import datetime as _dt
from typing import Any, Iterable

from .cache import HttpError, JsonFileCache, NetworkUnavailable
from .models import ChebiRelation, ChebiTerm, normalise_chebi_id

OLS4_TERM_URL = "https://www.ebi.ac.uk/ols4/api/v2/ontologies/chebi/classes/{iri}"
OLS4_ONTOLOGY_URL = "https://www.ebi.ac.uk/ols4/api/ontologies/chebi"
CHEBI_BACKEND_URL = "https://www.ebi.ac.uk/chebi/backend/api/public/compound/{num}/"
CACHE_NS = "chebi"
CACHE_NS_META = "chebi_meta"

OBO = "http://purl.obolibrary.org/obo/"
RO_HAS_ROLE = "RO_0000087"

# OLS4 property IRIs -> readable relation names
_RELATION_NAMES = {
    "RO_0000087": "has_role",
    "RO_0018036": "is_tautomer_of",
    "RO_0018033": "is_conjugate_base_of",
    "RO_0018034": "is_conjugate_acid_of",
    "BFO_0000051": "has_part",
    "chebi#has_functional_parent": "has_functional_parent",
    "chebi#has_parent_hydride": "has_parent_hydride",
    "chebi#is_enantiomer_of": "is_enantiomer_of",
    "chebi#is_substituent_group_from": "is_substituent_group_from",
}
INCOMING_KEPT = frozenset({"has part", "is conjugate acid of", "is conjugate base of", "is tautomer of", "has functional parent", "is enantiomer of", "has parent hydride"})
_BACKEND_RELATION_NAMES = {
    "is a": "is_a",
    "has role": "has_role",
    "is tautomer of": "is_tautomer_of",
    "is conjugate base of": "is_conjugate_base_of",
    "is conjugate acid of": "is_conjugate_acid_of",
    "has part": "has_part",
    "has functional parent": "has_functional_parent",
    "has parent hydride": "has_parent_hydride",
    "is enantiomer of": "is_enantiomer_of",
    "is substituent group from": "is_substituent_group_from",
}


def chebi_iri(chebi_id: str) -> str:
    return OBO + chebi_id.replace(":", "_")


def ols4_term_url(chebi_id: str) -> str:
    # OLS4 wants the IRI double URL-encoded inside the path
    from urllib.parse import quote

    return OLS4_TERM_URL.format(iri=quote(quote(chebi_iri(chebi_id), safe=""), safe=""))


def _iri_to_id(iri: str) -> str | None:
    tail = iri.rsplit("/", 1)[-1]
    return normalise_chebi_id(tail) if tail.startswith("CHEBI_") else None


def _first(v: Any) -> Any:
    if isinstance(v, list):
        return v[0] if v else None
    return v


def _now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat()


# ---------------------------------------------------------------------------
# Trimming of raw payloads (what gets cached / recorded as fixtures)
# ---------------------------------------------------------------------------

def trim_ols4(raw: dict[str, Any]) -> dict[str, Any]:
    keep = ["curie", "label", "definition", "directParent", "hierarchicalAncestor", "isObsolete", "numHierarchicalDescendants", "numDescendants", "https://w3id.org/chemrof/inchi_key_string", "http://www.geneontology.org/formats/oboInOwl#hasExactSynonym", "http://www.w3.org/2002/07/owl#deprecated", "http://purl.obolibrary.org/obo/IAO_0100001"]
    out = {k: raw[k] for k in keep if k in raw}
    out["relatedTo"] = [{"property": r.get("property"), "value": r.get("value")} for r in raw.get("relatedTo") or [] if isinstance(r, dict)]
    return out


def trim_backend(raw: dict[str, Any]) -> dict[str, Any]:
    out = {k: raw.get(k) for k in ("id", "chebi_accession", "name", "ascii_name", "definition", "stars", "secondary_ids", "is_released")}
    ds = raw.get("default_structure") or {}
    out["default_structure"] = {k: ds.get(k) for k in ("standard_inchi_key", "smiles")}
    rel = raw.get("ontology_relations") or {}
    out["ontology_relations"] = {
        "outgoing_relations": [{k: r.get(k) for k in ("relation_type", "final_id", "final_name")} for r in rel.get("outgoing_relations") or []],
        # incoming "is a"/"has role" edges are children/role-bearers and can number in the
        # thousands; only structural/family context (mixture has_part, conjugate forms, ...) is kept
        "incoming_relations": [{k: r.get(k) for k in ("relation_type", "init_id", "init_name")} for r in rel.get("incoming_relations") or [] if r.get("relation_type") in INCOMING_KEPT],
    }
    return out


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------

def _syn_values(v: Any) -> list[str]:
    out: list[str] = []
    for item in v or []:
        if isinstance(item, dict):
            val = item.get("value")
        else:
            val = item
        if isinstance(val, str) and val not in out:
            out.append(val)
    return out


def build_term(chebi_id: str, ols4: dict[str, Any] | None, backend: dict[str, Any] | None, release: str | None) -> ChebiTerm:
    term = ChebiTerm(chebi_id=chebi_id, chebi_release=release, retrieved_at=_now())
    if ols4:
        term.name = _first(ols4.get("label"))
        term.definition = _first(ols4.get("definition"))
        term.obsolete = bool(_first(ols4.get("isObsolete")) or _first(ols4.get("http://www.w3.org/2002/07/owl#deprecated")))
        term.is_a = sorted({i for i in (_iri_to_id(p) for p in ols4.get("directParent") or []) if i})
        term.ols4_ancestors = sorted({i for i in (_iri_to_id(p) for p in ols4.get("hierarchicalAncestor") or []) if i})
        nd = ols4.get("numHierarchicalDescendants")
        term.num_descendants = int(nd) if nd is not None else None
        term.inchikey = _first(ols4.get("https://w3id.org/chemrof/inchi_key_string"))
        term.synonyms = _syn_values(ols4.get("http://www.geneontology.org/formats/oboInOwl#hasExactSynonym")) + _syn_values(ols4.get("http://www.geneontology.org/formats/oboInOwl#hasRelatedSynonym"))
        roles: list[str] = []
        rels: list[ChebiRelation] = []
        for r in ols4.get("relatedTo") or []:
            prop = (r.get("property") or "").rsplit("/", 1)[-1]
            tgt = _iri_to_id(r.get("value") or "")
            if not tgt:
                continue
            if prop == RO_HAS_ROLE:
                if tgt not in roles:
                    roles.append(tgt)
            else:
                rels.append(ChebiRelation(relation=_RELATION_NAMES.get(prop, prop), target_id=tgt))
        term.roles = sorted(roles)
        term.relations = rels
        term.sources["ols4"] = ols4
    if backend:
        term.stars = backend.get("stars")
        if term.name is None:
            term.name = backend.get("ascii_name") or backend.get("name")
        if term.definition is None:
            term.definition = backend.get("definition")
        ds = backend.get("default_structure") or {}
        if not term.inchikey:
            term.inchikey = ds.get("standard_inchi_key")
        rel = backend.get("ontology_relations") or {}
        if not ols4:  # OLS4 missing the term: take graph edges from the backend
            term.is_a = sorted({f"CHEBI:{r['final_id']}" for r in rel.get("outgoing_relations") or [] if r.get("relation_type") == "is a"})
            term.roles = sorted({f"CHEBI:{r['final_id']}" for r in rel.get("outgoing_relations") or [] if r.get("relation_type") == "has role"})
            term.relations = [ChebiRelation(relation=_BACKEND_RELATION_NAMES.get(r["relation_type"], r["relation_type"]), target_id=f"CHEBI:{r['final_id']}", target_name=r.get("final_name")) for r in rel.get("outgoing_relations") or [] if r.get("relation_type") not in ("is a", "has role")]
        term.incoming = [ChebiRelation(relation=_BACKEND_RELATION_NAMES.get(r["relation_type"], r["relation_type"]), target_id=f"CHEBI:{r['init_id']}", target_name=r.get("init_name")) for r in rel.get("incoming_relations") or [] if r.get("init_id") is not None]
        term.sources["chebi_backend"] = backend
    return term


# ---------------------------------------------------------------------------
# Client / graph
# ---------------------------------------------------------------------------

class ChebiClient:
    def __init__(self, transport, cache: JsonFileCache, use_backend: bool = True):
        self.transport = transport
        self.cache = cache
        self.use_backend = use_backend
        self._terms: dict[str, ChebiTerm] = {}
        self._ancestors: dict[str, frozenset[str]] = {}
        self._release: str | None = None

    # -- release ---------------------------------------------------------------
    def release(self) -> str | None:
        if self._release is None:
            try:
                meta = self.cache.get_or_fetch(CACHE_NS_META, "ontology", lambda: {"version": self.transport.get_json(OLS4_ONTOLOGY_URL).get("version"), "recorded_at": _now()})
                self._release = str(meta.get("version")) if meta.get("version") is not None else "unknown"
            except (NetworkUnavailable, HttpError):
                self._release = "unknown"
        return self._release

    # -- raw fetches ------------------------------------------------------------
    def _fetch_record(self, chebi_id: str) -> dict[str, Any]:
        ols4: dict[str, Any] | None
        try:
            ols4 = trim_ols4(self.transport.get_json(ols4_term_url(chebi_id)))
        except HttpError as exc:
            if exc.status != 404:
                raise
            ols4 = None
        backend: dict[str, Any] | None = None
        if self.use_backend:
            try:
                backend = trim_backend(self.transport.get_json(CHEBI_BACKEND_URL.format(num=chebi_id.split(":")[1])))
            except HttpError as exc:
                if exc.status != 404:
                    raise
                backend = None
        if ols4 is None and backend is None:
            return {"chebi_id": chebi_id, "missing": True, "retrieved_at": _now()}
        return build_term(chebi_id, ols4, backend, self.release()).to_dict()

    # -- public ----------------------------------------------------------------
    def term(self, chebi_id: str) -> ChebiTerm:
        cid = normalise_chebi_id(chebi_id)
        if cid is None:
            raise ValueError(f"not a ChEBI id: {chebi_id!r}")
        if cid in self._terms:
            return self._terms[cid]
        rec = self.cache.get_or_fetch(CACHE_NS, cid, lambda: self._fetch_record(cid))
        if rec.get("missing"):
            t = ChebiTerm(chebi_id=cid, name=None, definition="(term not found in OLS4 or ChEBI backend)", retrieved_at=rec.get("retrieved_at"))
        else:
            t = ChebiTerm.from_dict(rec)
        self._terms[cid] = t
        return t

    def has_term(self, chebi_id: str) -> bool:
        t = self.term(chebi_id)
        return t.name is not None

    def name(self, chebi_id: str) -> str:
        return self.term(chebi_id).name or "?"

    def names(self, ids: Iterable[str]) -> dict[str, str]:
        return {i: self.name(i) for i in ids}

    def ancestors(self, chebi_id: str) -> frozenset[str]:
        """Transitive is_a closure, computed locally from direct parents (excludes the term itself)."""
        cid = normalise_chebi_id(chebi_id)
        if cid in self._ancestors:
            return self._ancestors[cid]
        seen: set[str] = set()
        stack = list(self.term(cid).is_a)
        while stack:
            p = stack.pop()
            if p in seen:
                continue
            seen.add(p)
            stack.extend(self.term(p).is_a)
        fs = frozenset(seen)
        self._ancestors[cid] = fs
        return fs

    def ancestor_depths(self, chebi_id: str) -> dict[str, int]:
        """Shortest is_a path length from ``chebi_id`` to each ancestor."""
        cid = normalise_chebi_id(chebi_id)
        depths: dict[str, int] = {}
        frontier = [(p, 1) for p in self.term(cid).is_a]
        while frontier:
            nxt = []
            for p, d in frontier:
                if p in depths and depths[p] <= d:
                    continue
                depths[p] = d
                nxt.extend((q, d + 1) for q in self.term(p).is_a)
            frontier = nxt
        return depths

    def direct_roles(self, chebi_id: str) -> frozenset[str]:
        return frozenset(self.term(chebi_id).roles)

    def inherited_roles(self, chebi_id: str) -> dict[str, list[str]]:
        """Roles asserted on the term itself or on any is_a ancestor -> the terms asserting them."""
        cid = normalise_chebi_id(chebi_id)
        out: dict[str, list[str]] = {}
        for t in [cid, *sorted(self.ancestors(cid))]:
            for r in self.term(t).roles:
                out.setdefault(r, []).append(t)
        return out

    def role_closure(self, role_ids: Iterable[str]) -> frozenset[str]:
        """Role ids plus all their is_a ancestors (e.g. antibacterial drug -> antibacterial agent -> antimicrobial agent)."""
        out: set[str] = set()
        for r in role_ids:
            out.add(r)
            out |= self.ancestors(r)
        return frozenset(out)

    def is_a_path(self, chebi_id: str, ancestor: str) -> list[str] | None:
        """One shortest is_a path from chebi_id up to ancestor (inclusive), or None."""
        cid = normalise_chebi_id(chebi_id)
        target = normalise_chebi_id(ancestor)
        prev: dict[str, str | None] = {cid: None}
        frontier = [cid]
        while frontier:
            nxt = []
            for t in frontier:
                if t == target:
                    path = []
                    cur: str | None = t
                    while cur is not None:
                        path.append(cur)
                        cur = prev[cur]
                    return list(reversed(path))
                for p in self.term(t).is_a:
                    if p not in prev:
                        prev[p] = t
                        nxt.append(p)
            frontier = nxt
        return None
