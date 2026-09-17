"""Wires the layers together: identity -> mapping -> ontology (-> classification, Phase 6)."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .cache import HttpTransport, JsonFileCache, NetworkUnavailable
from .chebi import ChebiClient
from .identity import IdentityNotFound, RcsbClient
import json

from .mapping import ChebiLookup, Mapper, UniChemClient, load_synonym_table
from .models import MAP_UNRESOLVED_NETWORK, MAP_UNRESOLVED_NOT_FOUND, ChemicalIdentity, MappingResult

DEFAULT_CACHE = Path("cache")
DEFAULT_SYNONYMS = Path("data/config/manual_synonyms.json")
DEFAULT_RELATED_PARENTS = Path("data/config/related_parents.json")

# ChEBI classes whose descendants are polymers/macromolecules/mixtures rather than a single
# molecular species; used only to rank equivalent ids sharing one InChIKey.
CLASS_LIKE_ANCESTORS = frozenset({"CHEBI:33839", "CHEBI:60027", "CHEBI:60004"})  # macromolecule, polymer, mixture

FAMILY_RELATIONS = ("has_part",)  # incoming: mixture/family  has_part  this component
FORM_RELATIONS = ("is_conjugate_acid_of", "is_conjugate_base_of", "is_tautomer_of")


@dataclass
class FamilyContext:
    """Mixture/family and ionisation-form relations of the exact deposited species."""

    member_of: list[dict[str, str]] = field(default_factory=list)   # e.g. gentamicin C1a -> gentamycin
    forms: list[dict[str, str]] = field(default_factory=list)       # conjugate acid/base, tautomers

    def to_dict(self) -> dict[str, Any]:
        return {"member_of": self.member_of, "forms": self.forms}


@dataclass
class Resolution:
    identity: ChemicalIdentity
    mapping: MappingResult
    family: FamilyContext = field(default_factory=FamilyContext)
    primary_name: str | None = None
    equivalent_names: dict[str, str] = field(default_factory=dict)


class Pipeline:
    def __init__(self, cache_dir: str | Path = DEFAULT_CACHE, transport=None, synonyms: str | Path | None = DEFAULT_SYNONYMS, related_parents: str | Path | None = DEFAULT_RELATED_PARENTS):
        self.cache = JsonFileCache(cache_dir)
        self.transport = transport or HttpTransport()
        self.rcsb = RcsbClient(self.transport, self.cache)
        self.unichem = UniChemClient(self.transport, self.cache)
        self.chebi = ChebiClient(self.transport, self.cache)
        self.mapper = Mapper(self.unichem, lookup=self.lookup, synonym_table=load_synonym_table(synonyms))
        self.related_parents: dict[str, dict[str, str]] = {}
        if related_parents and Path(related_parents).exists():
            self.related_parents = json.loads(Path(related_parents).read_text(encoding="utf-8"))

    def lookup(self, cid: str) -> ChebiLookup:
        """Facts about a ChEBI id used by the mapper. A monomer's InChIKey is reused by ChEBI for
        polymer classes (glucose -> glucans) and broad classes (glycerol -> alditol), so class-like
        entries are ranked last when choosing a primary id."""
        try:
            term = self.chebi.term(cid)
            anc = self.chebi.ancestors(cid)
        except NetworkUnavailable:
            return ChebiLookup()
        forms = {r.target_id for r in term.relations + term.incoming if r.relation in FORM_RELATIONS}
        return ChebiLookup(stars=term.stars, inchikey=term.inchikey, is_class_like=bool(anc & CLASS_LIKE_ANCESTORS), num_descendants=term.num_descendants, forms=frozenset(forms))

    def resolve(self, input_id: str) -> Resolution:
        try:
            ident = self.rcsb.identity(input_id)
        except IdentityNotFound:
            ident = ChemicalIdentity(input_id=input_id.strip().upper())
            return Resolution(ident, MappingResult(status=MAP_UNRESOLVED_NOT_FOUND, source_id=ident.input_id, notes=["id not found in the RCSB chemical component / BIRD dictionaries"]))
        except NetworkUnavailable as exc:
            ident = ChemicalIdentity(input_id=input_id.strip().upper())
            return Resolution(ident, MappingResult(status=MAP_UNRESOLVED_NETWORK, source_id=ident.input_id, notes=[f"identity retrieval failed: {exc}"]))
        mapping = self.mapper.map(ident)
        rp = self.related_parents.get(ident.input_id.upper())
        if rp:
            mapping.related_parent = dict(rp)
        res = Resolution(ident, mapping)
        if mapping.resolved:
            try:
                res.primary_name = self.chebi.name(mapping.primary_chebi_id)
                res.equivalent_names = self.chebi.names(mapping.equivalent_chebi_ids)
                res.family = self.family_context(mapping.primary_chebi_id)
            except NetworkUnavailable as exc:
                mapping.notes.append(f"ontology retrieval incomplete: {exc}")
        else:
            for c in mapping.candidates:
                try:
                    c.name = self.chebi.name(c.chebi_id)
                except NetworkUnavailable:
                    pass
        return res

    def family_context(self, chebi_id: str) -> FamilyContext:
        term = self.chebi.term(chebi_id)
        fc = FamilyContext()
        for r in term.incoming:
            if r.relation in FAMILY_RELATIONS:
                fc.member_of.append({"chebi_id": r.target_id, "name": r.target_name or self.chebi.name(r.target_id), "relation": f"{r.relation} (incoming)"})
        for r in term.relations + term.incoming:
            if r.relation in FORM_RELATIONS:
                fc.forms.append({"chebi_id": r.target_id, "name": r.target_name or self.chebi.name(r.target_id), "relation": r.relation})
        return fc


    # ------------------------------------------------------------------ ontology
    def ontology_summary(self, mapping: MappingResult) -> dict[str, Any]:
        """Complete ChEBI graph evidence for a resolved mapping (Phase 4 deliverable).

        Evidence is unioned over the equivalent ids (same standardised identity); each entry
        records which id it came from.
        """
        ids = mapping.chebi_ids_for_evidence
        anc: dict[str, dict[str, Any]] = {}
        direct_roles: dict[str, list[str]] = {}
        inherited: dict[str, dict[str, Any]] = {}
        for cid in ids:
            for a, d in self.chebi.ancestor_depths(cid).items():
                cur = anc.get(a)
                if cur is None or d < cur["min_depth"]:
                    anc[a] = {"chebi_id": a, "name": self.chebi.name(a), "min_depth": d, "via": cid}
            for r in self.chebi.direct_roles(cid):
                direct_roles.setdefault(r, []).append(cid)
            for r, via in self.chebi.inherited_roles(cid).items():
                inherited.setdefault(r, {"chebi_id": r, "name": self.chebi.name(r), "asserted_on": []})
                for v in via:
                    if v not in inherited[r]["asserted_on"]:
                        inherited[r]["asserted_on"].append(v)
        closure = self.chebi.role_closure(inherited.keys())
        implied = sorted(closure - set(inherited))
        return {
            "chebi_ids": ids,
            "terms": {cid: {"name": self.chebi.name(cid), "definition": self.chebi.term(cid).definition, "stars": self.chebi.term(cid).stars, "is_a": self.chebi.term(cid).is_a} for cid in ids},
            "is_a_ancestors": sorted(anc.values(), key=lambda x: (x["min_depth"], x["chebi_id"])),
            "direct_roles": {r: {"name": self.chebi.name(r), "on": v} for r, v in sorted(direct_roles.items())},
            "inherited_roles": dict(sorted(inherited.items())),
            "roles_implied_by_role_hierarchy": [{"chebi_id": r, "name": self.chebi.name(r)} for r in implied],
            "chebi_release": self.chebi.release(),
        }
