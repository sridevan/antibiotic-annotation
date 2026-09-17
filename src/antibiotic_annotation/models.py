"""Data models shared across identity retrieval, mapping, ontology and classification.

Everything is a plain dataclass with a ``to_dict`` so records can be written as JSON
and inspected without the package.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ENTITY_CCD = "ccd_nonpolymer"
ENTITY_BIRD = "bird_peptide_like"
ENTITY_CHEBI_INPUT = "chebi_direct"
ENTITY_UNKNOWN = "unknown"

# Mapping statuses
MAP_RESOLVED = "resolved"
MAP_UNRESOLVED_NO_CHEBI = "unresolved_no_chebi"
MAP_UNRESOLVED_CONFLICT = "unresolved_identity_conflict"
MAP_UNRESOLVED_NO_STRUCTURE = "unresolved_no_structure"
MAP_UNRESOLVED_NETWORK = "unresolved_network"
MAP_UNRESOLVED_NOT_FOUND = "unresolved_identity_not_found"

# Mapping methods
METHOD_DIRECT_XREF = "direct_ccd_crossref"
METHOD_UNICHEM_INCHIKEY = "unichem_inchikey"
METHOD_UNICHEM_CONNECTIVITY = "unichem_connectivity_protonation"
METHOD_MANUAL_SYNONYM = "manual_synonym"
METHOD_CHEBI_INPUT = "chebi_input"
METHOD_UNRESOLVED = "unresolved"


def _clean(d: Any) -> Any:
    """Recursively convert dataclasses/sets to JSON-friendly structures."""
    if hasattr(d, "to_dict"):
        return d.to_dict()
    if isinstance(d, dict):
        return {k: _clean(v) for k, v in d.items()}
    if isinstance(d, (list, tuple)):
        return [_clean(v) for v in d]
    if isinstance(d, (set, frozenset)):
        return sorted(_clean(v) for v in d)
    return d


@dataclass
class ChemicalIdentity:
    """What the PDB (RCSB) says the deposited chemical species is."""

    input_id: str
    entity_kind: str = ENTITY_UNKNOWN
    name: str | None = None
    ccd_type: str | None = None
    formula: str | None = None
    formal_charge: int | None = None
    inchi: str | None = None
    inchikey: str | None = None
    smiles: str | None = None
    synonyms: list[str] = field(default_factory=list)
    xrefs: dict[str, list[str]] = field(default_factory=dict)
    atc_codes: list[str] = field(default_factory=list)
    bird_class: str | None = None
    bird_type: str | None = None
    source: str = "RCSB"
    retrieved_at: str | None = None

    @property
    def inchikey_connectivity(self) -> str | None:
        return self.inchikey.split("-")[0] if self.inchikey else None

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["inchikey_connectivity"] = self.inchikey_connectivity
        return d


@dataclass
class MappingCandidate:
    """A ChEBI entity considered during mapping, with how it compares to the query."""

    chebi_id: str
    inchikey: str | None = None
    name: str | None = None
    comparison: dict[str, bool] = field(default_factory=dict)
    mismatches: list[str] = field(default_factory=list)
    accepted: bool = False
    note: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class MappingResult:
    """How a chemical identity was mapped to ChEBI, and how confident that is."""

    status: str = MAP_UNRESOLVED_NO_CHEBI
    method: str = METHOD_UNRESOLVED
    source_id: str | None = None  # identifier used for the lookup (InChIKey, CCD id, ...)
    primary_chebi_id: str | None = None
    equivalent_chebi_ids: list[str] = field(default_factory=list)  # same standard InChIKey
    evidence_unioned: bool = False
    candidates: list[MappingCandidate] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def resolved(self) -> bool:
        return self.status == MAP_RESOLVED and self.primary_chebi_id is not None

    @property
    def chebi_ids_for_evidence(self) -> list[str]:
        if not self.resolved:
            return []
        ids = [self.primary_chebi_id] + [c for c in self.equivalent_chebi_ids if c != self.primary_chebi_id]
        return ids

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["resolved"] = self.resolved
        return d


@dataclass
class ChebiRelation:
    relation: str  # e.g. "is_a", "has_role", "is_tautomer_of", "has_part"
    target_id: str
    target_name: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class ChebiTerm:
    """Normalised ChEBI term merged from OLS4 (graph) and the ChEBI backend (curation data)."""

    chebi_id: str
    name: str | None = None
    definition: str | None = None
    stars: int | None = None
    obsolete: bool = False
    is_a: list[str] = field(default_factory=list)  # direct parents
    roles: list[str] = field(default_factory=list)  # direct has_role targets
    relations: list[ChebiRelation] = field(default_factory=list)  # other outgoing relations
    incoming: list[ChebiRelation] = field(default_factory=list)  # relations pointing at this term
    synonyms: list[str] = field(default_factory=list)
    inchikey: str | None = None
    ols4_ancestors: list[str] = field(default_factory=list)  # OLS4 precomputed, cross-check only
    num_descendants: int | None = None
    chebi_release: str | None = None
    retrieved_at: str | None = None
    sources: dict[str, Any] = field(default_factory=dict)  # trimmed raw payloads

    def to_dict(self) -> dict[str, Any]:
        return _clean(asdict(self))

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "ChebiTerm":
        d = dict(d)
        d["relations"] = [ChebiRelation(**r) for r in d.get("relations", [])]
        d["incoming"] = [ChebiRelation(**r) for r in d.get("incoming", [])]
        return cls(**d)


@dataclass
class BenchmarkItem:
    compound: str
    label: int | None  # 1 antibiotic, 0 not, None = no truth (challenge / peptide track)
    benchmark_set: str  # primary_positive | primary_negative | paper_positive | challenge | peptide_track
    stratum: str
    input_id: str | None = None  # CCD, PRD or CHEBI id used for annotation
    entity_kind: str | None = None
    pdb_examples: list[str] = field(default_factory=list)
    workbook_chebi_id: str | None = None
    note: str | None = None
    source_sheet: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def normalise_chebi_id(value: str | int | None) -> str | None:
    """Return 'CHEBI:12345' for any of 12345, 'CHEBI_12345', 'chebi:12345', 'CHEBI:12345'."""
    if value is None:
        return None
    s = str(value).strip()
    if not s:
        return None
    s = s.replace("CHEBI_", "CHEBI:").replace("chebi:", "CHEBI:").replace("Chebi:", "CHEBI:")
    if s.isdigit():
        return f"CHEBI:{s}"
    if s.upper().startswith("CHEBI:"):
        num = s.split(":", 1)[1]
        if num.isdigit():
            return f"CHEBI:{num}"
    return None
