"""Wires the layers together: identity -> mapping -> ontology (-> classification, Phase 6)."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .cache import HttpTransport, JsonFileCache, NetworkUnavailable
from .chebi import ChebiClient
from .identity import IdentityNotFound, RcsbClient
from .mapping import Mapper, UniChemClient, load_synonym_table
from .models import MAP_UNRESOLVED_NETWORK, MAP_UNRESOLVED_NOT_FOUND, ChemicalIdentity, MappingResult

DEFAULT_CACHE = Path("cache")
DEFAULT_SYNONYMS = Path("data/config/manual_synonyms.json")

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
    def __init__(self, cache_dir: str | Path = DEFAULT_CACHE, transport=None, synonyms: str | Path | None = DEFAULT_SYNONYMS):
        self.cache = JsonFileCache(cache_dir)
        self.transport = transport or HttpTransport()
        self.rcsb = RcsbClient(self.transport, self.cache)
        self.unichem = UniChemClient(self.transport, self.cache)
        self.chebi = ChebiClient(self.transport, self.cache)
        self.mapper = Mapper(
            self.unichem,
            rank=self.rank_key,
            inchikey_of=lambda cid: self._safe_inchikey(cid),
            synonym_table=load_synonym_table(synonyms),
        )

    def rank_key(self, cid: str) -> tuple:
        """Sort key for choosing the primary among ChEBI ids that share a standard InChIKey.

        A monomer's InChIKey is reused by ChEBI for polymer classes (glucose -> glucans) and
        broad classes (glycerol -> alditol), so prefer: not a macromolecule/mixture/polymer
        class, fewest is_a descendants, highest star rating. Smaller tuple = preferred.
        """
        try:
            term = self.chebi.term(cid)
            anc = self.chebi.ancestors(cid)
        except NetworkUnavailable:
            return (1, 10**6, 0)
        is_class = int(bool(anc & CLASS_LIKE_ANCESTORS))
        return (is_class, term.num_descendants if term.num_descendants is not None else 10**6, -(term.stars or 0))

    def _safe_inchikey(self, cid: str) -> str | None:
        try:
            return self.chebi.term(cid).inchikey
        except NetworkUnavailable:
            return None

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
