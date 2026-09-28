"""Wires the layers together: compound identity -> ChEBI mapping -> ontology context."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .cache import HttpTransport, JsonFileCache, NetworkUnavailable
from .chebi import ChebiClient
from .identity import CompoundClient, IdentityNotFound
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
        self.compounds = CompoundClient(self.transport, self.cache)
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
            ident = self.compounds.identity(input_id)
        except (IdentityNotFound, NetworkUnavailable) as exc:
            status, note = (MAP_UNRESOLVED_NOT_FOUND, "id not found in the PDBe compound (CCD / BIRD) dictionaries") if isinstance(exc, IdentityNotFound) else (MAP_UNRESOLVED_NETWORK, f"identity retrieval failed: {exc}")
            ident = ChemicalIdentity(input_id=input_id.strip().upper())
            return Resolution(ident, MappingResult(status=status, source_id=ident.input_id, notes=[note]))
        mapping = self.mapper.map(ident)
        mapping.related_parent = self.related_parents.get(ident.input_id.upper())
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

        def entry(r, relation: str) -> dict[str, str]:
            return {"chebi_id": r.target_id, "name": r.target_name or self.chebi.name(r.target_id), "relation": relation}

        forms = {(r.target_id, r.relation): r for r in term.relations + term.incoming if r.relation in FORM_RELATIONS}  # dedupe
        return FamilyContext(
            member_of=[entry(r, f"{r.relation} (incoming)") for r in term.incoming if r.relation in FAMILY_RELATIONS],
            forms=[entry(r, r.relation) for r in forms.values()],
        )
