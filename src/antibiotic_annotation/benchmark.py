"""Benchmark construction from the curated workbook, used for validation only (the production API never reads it).

The workbook is the source of truth for compound names, labels and strata. The mapping from
compound name to the deposited PDB chemical species (CCD id, PRD id, or a ChEBI id for
controls without a single CCD) lives in ``data/benchmark/benchmark_ids.tsv``; it was built by
resolving the workbook's "PDB examples" column through the RCSB entry API, not by name
matching, and is checked in so the benchmark is reproducible and auditable.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import openpyxl

from .models import ENTITY_BIRD, ENTITY_CCD, ENTITY_CHEBI_INPUT, BenchmarkItem, normalise_chebi_id

SET_PRIMARY_POS = "primary_positive"
SET_PRIMARY_NEG = "primary_negative"
SET_PAPER_POS = "paper_positive"
SET_CHALLENGE = "challenge"
SET_PEPTIDE = "peptide_track"

SHEET_MANUAL = "Manual benchmark audit"
SHEET_NEGATIVE = "Negative controls"
SHEET_PAPER = "Paper benchmark"

ENTITY_KIND_POLYMER_PEPTIDE = "polymer_peptide"


@dataclass
class IdTableRow:
    compound: str
    input_id: str
    entity_kind: str
    stratum: str
    resolved_from_pdb_entry: str
    note: str


@dataclass
class BenchmarkSet:
    primary_positives: list[BenchmarkItem] = field(default_factory=list)
    primary_negatives: list[BenchmarkItem] = field(default_factory=list)
    paper_positives: list[BenchmarkItem] = field(default_factory=list)
    challenge: list[BenchmarkItem] = field(default_factory=list)
    peptide_track: list[BenchmarkItem] = field(default_factory=list)

    @property
    def primary_items(self) -> list[BenchmarkItem]:
        return self.primary_positives + self.primary_negatives

    @property
    def all_items(self) -> list[BenchmarkItem]:
        return self.primary_items + self.paper_positives + self.challenge + self.peptide_track

    @property
    def annotatable_items(self) -> list[BenchmarkItem]:
        """Everything that goes through the small-molecule/BIRD pipeline (peptide track excluded)."""
        return self.primary_items + self.paper_positives + self.challenge


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rows(ws) -> list[dict[str, str | None]]:
    rows = list(ws.iter_rows(values_only=True))
    header = [str(h).strip() if h is not None else "" for h in rows[0]]
    out = []
    for r in rows[1:]:
        if not any(c is not None for c in r):
            continue
        d = {header[i]: (str(r[i]).strip() if i < len(r) and r[i] is not None else None) for i in range(len(header))}
        out.append(d)
    return out


def _norm_name(s: str) -> str:
    return " ".join(s.lower().replace("-", " ").split())


def load_id_table(path: str | Path) -> dict[str, IdTableRow]:
    table: dict[str, IdTableRow] = {}
    with Path(path).open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            row = IdTableRow(**{k: (v or "") for k, v in r.items()})
            table[_norm_name(row.compound)] = row
    return table


def _split_pdb_examples(s: str | None) -> list[str]:
    if not s:
        return []
    return [p.strip().upper() for p in s.split(",") if p.strip()]


def _negative_stratum(control_stratum: str | None, expected_label: str | None) -> tuple[str, int | None, str]:
    """Return (stratum, label, benchmark_set) for a Negative-controls row."""
    cs = (control_stratum or "").lower()
    if cs.startswith("easy negative"):
        return "easy", 0, SET_PRIMARY_NEG
    if cs.startswith("hard negative"):
        return "hard_ontology", 0, SET_PRIMARY_NEG
    if cs.startswith("challenge negative"):
        return "challenge_negative", 0, SET_PRIMARY_NEG
    if cs.startswith("challenge / policy") or "flag for review" in (expected_label or "").lower():
        return "challenge_policy", None, SET_CHALLENGE
    if cs.startswith("challenge:"):
        # e.g. the second fluconazole row: a clean negative for this project's definition
        return "hard_ontology", 0, SET_PRIMARY_NEG
    raise ValueError(f"unrecognised control stratum: {control_stratum!r}")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def load_benchmark(workbook: str | Path, id_table: str | Path | None = None) -> BenchmarkSet:
    workbook = Path(workbook)
    if id_table is None:
        id_table = workbook.parent / "benchmark_ids.tsv"
    ids = load_id_table(id_table)
    wb = openpyxl.load_workbook(workbook, read_only=True, data_only=True)
    bench = BenchmarkSet()
    seen_ids: set[str] = set()

    # ---- Manual benchmark audit -------------------------------------------------
    for r in _rows(wb[SHEET_MANUAL]):
        compound = r["Compound"]
        if not compound:
            continue
        key = _norm_name(compound)
        if key not in ids:
            raise KeyError(f"{compound!r} from sheet {SHEET_MANUAL!r} has no row in the id table")
        row = ids[key]
        entity_type = (r.get("Entity type") or "").lower()
        audit = (r.get("Audit outcome") or "").lower()
        item = BenchmarkItem(
            compound=compound,
            label=None,
            benchmark_set="",
            stratum=row.stratum,
            input_id=row.input_id,
            entity_kind=row.entity_kind,
            pdb_examples=_split_pdb_examples(r.get("PDB examples")),
            workbook_chebi_id=normalise_chebi_id(r.get("ChEBI ID")),
            note=r.get("Why it matters"),
            source_sheet=SHEET_MANUAL,
        )
        if row.entity_kind == ENTITY_KIND_POLYMER_PEPTIDE:
            item.benchmark_set = SET_PEPTIDE
            # Regulatory peptides are explicitly "Not antibiotic"; the AMP has no truth in this track.
            item.label = 0 if "not antibiotic" in audit else None
            bench.peptide_track.append(item)
        elif "not antibiotic" in audit or "metabolite" in entity_type:
            item.benchmark_set = SET_PRIMARY_NEG
            item.label = 0
            item.stratum = "same_domain"
            # Ornithine is also in the Negative-controls sheet under the same CCD; keep one row.
            if item.input_id not in seen_ids:
                bench.primary_negatives.append(item)
                seen_ids.add(item.input_id)
        else:
            item.benchmark_set = SET_PRIMARY_POS
            item.label = 1
            bench.primary_positives.append(item)
            seen_ids.add(item.input_id)

    # ---- Negative controls ------------------------------------------------------
    for r in _rows(wb[SHEET_NEGATIVE]):
        control = r["Control"]
        if not control:
            continue
        stratum, label, bset = _negative_stratum(r.get("Control stratum"), r.get("Expected label"))
        ccd = (r.get("PDB CCD") or "").strip()
        key = _norm_name(control.split("(")[0])
        if key in ids:
            row = ids[key]
            input_id, kind = row.input_id, row.entity_kind
        elif ccd and " " not in ccd:
            input_id, kind = ccd.upper(), ENTITY_CCD
        else:
            raise KeyError(f"cannot determine an input id for control {control!r}")
        if input_id in seen_ids:
            # duplicate deposited species (second fluconazole row, ornithine): keep the first
            continue
        seen_ids.add(input_id)
        item = BenchmarkItem(
            compound=control.split("(")[0].strip(),
            label=label,
            benchmark_set=bset,
            stratum=stratum,
            input_id=input_id,
            entity_kind=kind,
            workbook_chebi_id=normalise_chebi_id(r.get("ChEBI ID")),
            note=r.get("Ontology trap / signal"),
            source_sheet=SHEET_NEGATIVE,
        )
        # Same-domain ribosome ligands are still easy negatives for the split; keep sheet stratum.
        if bset == SET_PRIMARY_NEG:
            bench.primary_negatives.append(item)
        else:
            bench.challenge.append(item)

    # ---- Paper benchmark (external validation) ---------------------------------
    primary_ids = {i.input_id for i in bench.primary_positives}
    for r in _rows(wb[SHEET_PAPER]):
        compound = r["Compound"]
        ccd = (r.get("PDB CCD code reported in paper/tables") or "").strip().upper()
        if not compound or not ccd:
            continue
        if ccd in primary_ids:
            continue  # same deposited species already in the primary set
        key = _norm_name(compound)
        row = ids.get(key)
        stratum = row.stratum if row else _paper_class_to_stratum(r.get("Class in paper"))
        input_id, kind, note = ccd, ENTITY_CCD, f"Class in paper: {r.get('Class in paper')}"
        if row and row.resolved_from_pdb_entry.startswith("paper") and row.input_id.upper() != ccd:
            # audited correction of a wrong identifier in the workbook (see benchmark_ids.tsv note);
            # only id-table rows sourced from the paper sheet may override the paper's CCD code
            input_id, kind = row.input_id, row.entity_kind
            note += f" | benchmark identifier corrected from {ccd} to {row.input_id}: {row.note}"
            if input_id in primary_ids:
                continue
        item = BenchmarkItem(
            compound=compound,
            label=1,
            benchmark_set=SET_PAPER_POS,
            stratum=stratum,
            input_id=input_id,
            entity_kind=kind,
            note=note,
            source_sheet=SHEET_PAPER,
        )
        bench.paper_positives.append(item)

    return bench


def _paper_class_to_stratum(cls: str | None) -> str:
    c = (cls or "").lower()
    if "tetracycl" in c:
        return "tetracycline"
    if "aminoglycoside" in c or "aminocyclitol" in c:
        return "aminoglycoside"
    if "tuberactinomycin" in c:
        return "peptide_like"
    if "orthosomycin" in c:
        return "orthosomycin"
    if "pleuromutilin" in c:
        return "pleuromutilin"
    if "lincosamide" in c:
        return "lincosamide"
    return "other"


# ---------------------------------------------------------------------------
# Mapping run over the benchmark (Phase 3 report)
# ---------------------------------------------------------------------------

MAPPING_REPORT_COLUMNS = [
    "compound", "benchmark_set", "label", "stratum", "input_id", "entity_kind", "deposited_name", "formal_charge",
    "inchikey", "mapping_status", "mapping_method", "mapping_confidence", "primary_chebi_id", "primary_name", "equivalent_chebi_ids",
    "evidence_unioned", "family_context", "related_parent", "candidates", "workbook_chebi_id", "workbook_agreement", "bird_class", "notes",
]


def run_mapping(bench: "BenchmarkSet", pipeline) -> list[dict]:
    rows: list[dict] = []
    for item in bench.annotatable_items:
        res = pipeline.resolve(item.input_id)
        ident, m = res.identity, res.mapping
        if m.resolved:
            if item.workbook_chebi_id is None:
                agreement = "workbook had no id"
            elif item.workbook_chebi_id in m.equivalent_chebi_ids:
                agreement = "agrees"
            else:
                agreement = f"differs (workbook {item.workbook_chebi_id})"
        else:
            agreement = "unresolved" if item.workbook_chebi_id is None else f"unresolved (workbook {item.workbook_chebi_id})"
        rows.append({
            "compound": item.compound,
            "benchmark_set": item.benchmark_set,
            "label": "" if item.label is None else item.label,
            "stratum": item.stratum,
            "input_id": item.input_id,
            "entity_kind": ident.entity_kind,
            "deposited_name": ident.name or "",
            "formal_charge": "" if ident.formal_charge is None else ident.formal_charge,
            "inchikey": ident.inchikey or "",
            "mapping_status": m.status,
            "mapping_method": m.method,
            "mapping_confidence": m.confidence,
            "primary_chebi_id": m.primary_chebi_id or "",
            "primary_name": res.primary_name or "",
            "equivalent_chebi_ids": ";".join(f"{c} ({res.equivalent_names.get(c, '?')})" for c in m.equivalent_chebi_ids),
            "evidence_unioned": m.evidence_unioned,
            "family_context": ";".join(f"{f['name']} ({f['chebi_id']})" for f in res.family.member_of),
            "related_parent": f"{m.related_parent.get('name')} ({m.related_parent.get('chebi_id')}): {m.related_parent.get('relation')}" if m.related_parent else "",
            "candidates": ";".join(f"{c.chebi_id} {c.name or ''} [{'accepted' if c.accepted else 'rejected'}: {','.join(c.mismatches) or 'exact'}]" for c in m.candidates),
            "workbook_chebi_id": item.workbook_chebi_id or "",
            "workbook_agreement": agreement,
            "bird_class": ident.bird_class or "",
            "notes": " | ".join(m.notes),
            "_resolution": res,
        })
    return rows


def write_tsv(rows: list[dict], path: str | Path, columns: list[str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
