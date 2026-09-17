"""Public API.

    find_antibiotic_entities("5J7L", 1)          -> list[AntibioticEntity]
    inspect_assembly("5J7L", 1)                  -> AssemblyResult (hits + diagnostics)
    annotate_entity("TAC", "CCD")                -> AnnotationRecord
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .assembly import KIND_CCD, KIND_PRD, AssemblyClient, AssemblyEntity
from .classifier import Evidence, chebi_evidence, decide, name_stem_flags
from .models import ENTITY_BIRD, ENTITY_CCD
from .pipeline import Pipeline

STATUS_HIT = "antibiotic_like"
STATUS_NOT = "not_antibiotic"


@dataclass
class AnnotationRecord:
    entity_id: str
    entity_kind: str                     # "CCD" | "PRD"
    name: str | None
    chebi_id: str | None
    chebi_name: str | None
    antibiotic_like: bool
    evidence: Evidence
    reason: list[str]
    mapping: dict[str, Any]
    status: str                          # antibiotic_like | not_antibiotic | <mapping status when unresolved>
    family: dict[str, Any] = field(default_factory=dict)
    name_flags: dict[str, Any] = field(default_factory=dict)  # diagnostic only: naming stems such as -mycin

    def to_dict(self) -> dict[str, Any]:
        return {
            "entity_id": self.entity_id, "entity_kind": self.entity_kind, "name": self.name,
            "chebi_id": self.chebi_id, "chebi_name": self.chebi_name, "antibiotic_like": self.antibiotic_like,
            "evidence": self.evidence.to_dict(), "reason": self.reason, "mapping": self.mapping, "status": self.status, "family": self.family,
            "name_flags": self.name_flags,
        }


@dataclass
class AntibioticEntity(AnnotationRecord):
    pdb_id: str = ""
    assembly_id: str = ""
    chains: list[str] = field(default_factory=list)
    copies: int = 0

    def to_dict(self) -> dict[str, Any]:
        d = {"pdb_id": self.pdb_id, "assembly_id": self.assembly_id}
        d.update(super().to_dict())
        d["chains"] = self.chains
        d["copies"] = self.copies
        return d


@dataclass
class AssemblyResult:
    pdb_id: str
    assembly_id: str
    entities_inspected: list[AssemblyEntity]
    hits: list[AntibioticEntity]
    diagnostics: list[dict[str, Any]]     # every inspected entity with its status (hits included)

    def to_dict(self) -> dict[str, Any]:
        return {"pdb_id": self.pdb_id, "assembly_id": self.assembly_id, "entities_inspected": [e.to_dict() for e in self.entities_inspected], "hits": [h.to_dict() for h in self.hits], "diagnostics": self.diagnostics}


_default_pipeline: Pipeline | None = None


def _pipe(pipeline: Pipeline | None) -> Pipeline:
    global _default_pipeline
    if pipeline is not None:
        return pipeline
    if _default_pipeline is None:
        _default_pipeline = Pipeline()
    return _default_pipeline


def _input_id(entity_id: str, entity_kind: str) -> str:
    kind = entity_kind.upper()
    if kind not in (KIND_CCD, KIND_PRD):
        raise ValueError(f"entity_kind must be 'CCD' or 'PRD', got {entity_kind!r}")
    eid = entity_id.strip().upper()
    if kind == KIND_PRD and not eid.startswith("PRD_"):
        raise ValueError(f"PRD entity id must look like PRD_000226, got {entity_id!r}")
    if kind == KIND_CCD and eid.startswith("PRD_"):
        raise ValueError(f"{entity_id!r} is a PRD id, not a CCD code")
    return eid


def annotate_entity(entity_id: str, entity_kind: str, pipeline: Pipeline | None = None) -> AnnotationRecord:
    pipe = _pipe(pipeline)
    kind = entity_kind.upper()
    res = pipe.resolve(_input_id(entity_id, kind))
    ident, m = res.identity, res.mapping
    ev = Evidence()
    if ident.entity_kind == ENTITY_BIRD:
        ev.bird_class = ident.bird_class
        ev.bird_antibiotic = (ident.bird_class or "").strip().lower() == "antibiotic"
    if m.resolved:
        chebi_ev = chebi_evidence(pipe.chebi, m.chebi_ids_for_evidence)
        chebi_ev.bird_class, chebi_ev.bird_antibiotic = ev.bird_class, ev.bird_antibiotic
        ev = chebi_ev
    hit, reasons = decide(ev)
    if hit:
        status = STATUS_HIT
    elif m.resolved:
        status = STATUS_NOT
    else:
        status = m.status
    mapping = {
        "status": m.status, "method": m.method, "confidence": m.confidence,
        "equivalent_chebi_ids": m.equivalent_chebi_ids, "evidence_unioned": m.evidence_unioned,
        "candidates": [c.to_dict() for c in m.candidates], "related_parent": m.related_parent, "notes": m.notes,
        "inchikey": ident.inchikey,
    }
    names: dict[str, list[str]] = {"deposited_name": [ident.name] if ident.name else [], "deposited_synonyms": ident.synonyms}
    if m.resolved:
        names["chebi_name"] = [res.primary_name] if res.primary_name else []
        names["chebi_synonyms"] = [syn for cid in m.chebi_ids_for_evidence for syn in pipe.chebi.term(cid).synonyms]
    return AnnotationRecord(
        entity_id=ident.input_id, entity_kind=kind, name=res.primary_name or ident.name,
        chebi_id=m.primary_chebi_id, chebi_name=res.primary_name, antibiotic_like=hit, evidence=ev,
        reason=reasons, mapping=mapping, status=status, family=res.family.to_dict(), name_flags=name_stem_flags(names),
    )


def inspect_assembly(pdb_id: str, assembly_id: str | int, pipeline: Pipeline | None = None) -> AssemblyResult:
    """Annotate every CCD/PRD entity of one biological assembly; hits plus per-entity diagnostics."""
    pipe = _pipe(pipeline)
    client = AssemblyClient(pipe.transport, pipe.cache)
    pid, aid = pdb_id.strip().upper(), str(assembly_id).strip()
    entities = client.entities(pid, aid)
    hits: list[AntibioticEntity] = []
    diagnostics: list[dict[str, Any]] = []
    for ent in entities:
        rec = annotate_entity(ent.entity_id, ent.entity_kind, pipeline=pipe)
        diag = {"entity_id": ent.entity_id, "entity_kind": ent.entity_kind, "name": rec.name, "status": rec.status, "chebi_id": rec.chebi_id, "mapping_method": rec.mapping["method"],
                "antibiotic_naming_stem": rec.name_flags["antibiotic_naming_stem"], "naming_stems": sorted({x["stem"] for x in rec.name_flags["matches"]})}
        if rec.reason:
            diag["reason"] = rec.reason
        diagnostics.append(diag)
        if rec.antibiotic_like:
            hits.append(AntibioticEntity(**rec.__dict__, pdb_id=pid, assembly_id=aid, chains=ent.chains, copies=ent.copies))
    return AssemblyResult(pdb_id=pid, assembly_id=aid, entities_inspected=entities, hits=hits, diagnostics=diagnostics)


def find_antibiotic_entities(pdb_id: str, assembly_id: str | int, pipeline: Pipeline | None = None) -> list[AntibioticEntity]:
    return inspect_assembly(pdb_id, assembly_id, pipeline=pipeline).hits
