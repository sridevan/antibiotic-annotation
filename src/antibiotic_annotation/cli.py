"""Command-line interface.

    python -m antibiotic_annotation.cli annotate TAC
    python -m antibiotic_annotation.cli map-benchmark --input data/benchmark/antibiotic_ontology_rule_evaluation.xlsx
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

from .benchmark import MAPPING_REPORT_COLUMNS, load_benchmark, run_mapping, write_tsv
from .pipeline import Pipeline


def cmd_annotate(args: argparse.Namespace) -> int:
    pipe = Pipeline(cache_dir=args.cache)
    res = pipe.resolve(args.id)
    ident, m = res.identity, res.mapping
    if args.json:
        print(json.dumps({"identity": ident.to_dict(), "mapping": m.to_dict(), "primary_name": res.primary_name, "family": res.family.to_dict()}, indent=2))
        return 0
    print(ident.input_id)
    print(ident.name or "(no name)")
    print(f"entity_kind: {ident.entity_kind}")
    if ident.inchikey:
        print(f"inchikey: {ident.inchikey}")
    print(f"mapping: {m.status} via {m.method}")
    if m.resolved:
        print(f"chebi: {m.primary_chebi_id} ({res.primary_name})")
        if len(m.equivalent_chebi_ids) > 1:
            print("equivalent ids (same standardised identity, evidence unioned): " + ", ".join(f"{c} ({res.equivalent_names.get(c, '?')})" for c in m.equivalent_chebi_ids))
        for f in res.family.member_of:
            print(f"family: component of {f['name']} ({f['chebi_id']})")
    for c in m.candidates:
        if not c.accepted:
            print(f"candidate not accepted: {c.chebi_id} {c.name or ''} differs in {','.join(c.mismatches)}")
    for n in m.notes:
        print(f"note: {n}")
    return 0


def cmd_map_benchmark(args: argparse.Namespace) -> int:
    bench = load_benchmark(args.input)
    pipe = Pipeline(cache_dir=args.cache)
    rows = run_mapping(bench, pipe)
    out = Path(args.output)
    write_tsv(rows, out / "mapping_report.tsv", MAPPING_REPORT_COLUMNS)
    (out / "mapping").mkdir(parents=True, exist_ok=True)
    for r in rows:
        res = r["_resolution"]
        rec = {"identity": res.identity.to_dict(), "mapping": res.mapping.to_dict(), "primary_name": res.primary_name, "equivalent_names": res.equivalent_names, "family": res.family.to_dict()}
        (out / "mapping" / f"{r['input_id'].replace(':', '_')}.json").write_text(json.dumps(rec, indent=2, ensure_ascii=False))
    by_set = Counter((r["benchmark_set"], r["mapping_status"]) for r in rows)
    print(f"wrote {out / 'mapping_report.tsv'} ({len(rows)} compounds)")
    for (bset, status), n in sorted(by_set.items()):
        print(f"  {bset:18} {status:32} {n}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="antibiotic_annotation")
    p.add_argument("--cache", default="cache", help="cache directory (default: cache/)")
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("annotate", help="annotate one CCD / PRD / ChEBI id")
    a.add_argument("id")
    a.add_argument("--json", action="store_true")
    a.set_defaults(func=cmd_annotate)
    b = sub.add_parser("map-benchmark", help="map every benchmark compound and write output/mapping_report.tsv")
    b.add_argument("--input", default="data/benchmark/antibiotic_ontology_rule_evaluation.xlsx")
    b.add_argument("--output", default="output")
    b.set_defaults(func=cmd_map_benchmark)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
