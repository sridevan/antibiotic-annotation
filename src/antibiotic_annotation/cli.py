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

from .benchmark import MAPPING_REPORT_COLUMNS, load_benchmark, run_mapping, stratified_split, write_tsv
from .coverage import COVERAGE_COLUMNS, build_term_coverage
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


def _load_or_make_split(bench, out: Path, seed: str) -> dict[str, str]:
    path = out / "split.json"
    if path.exists():
        return json.loads(path.read_text())["assignment"]
    assignment = stratified_split(bench.primary_items, seed=seed, dev_fraction=0.7)
    out.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"seed": seed, "dev_fraction": 0.7, "assignment": assignment, "items": {i.input_id: {"compound": i.compound, "label": i.label, "stratum": i.stratum} for i in bench.primary_items}}, indent=1, sort_keys=True))
    return assignment


def cmd_split(args: argparse.Namespace) -> int:
    bench = load_benchmark(args.input)
    out = Path(args.output)
    if (out / "split.json").exists() and not args.force:
        print(f"{out / 'split.json'} already exists (use --force to regenerate)")
        return 0
    if args.force and (out / "split.json").exists():
        (out / "split.json").unlink()
    assignment = _load_or_make_split(bench, out, args.seed)
    for part in ("dev", "holdout"):
        items = [i for i in bench.primary_items if assignment[i.input_id] == part]
        print(f"{part}: {len(items)} items, {sum(i.label == 1 for i in items)} positives, {sum(i.label == 0 for i in items)} negatives")
        for i in sorted(items, key=lambda x: (-(x.label or 0), x.stratum, x.compound)):
            print(f"   {i.label} {i.stratum:18} {i.compound} ({i.input_id})")
    return 0


def cmd_coverage(args: argparse.Namespace) -> int:
    bench = load_benchmark(args.input)
    out = Path(args.output)
    assignment = _load_or_make_split(bench, out, args.seed)
    pipe = Pipeline(cache_dir=args.cache)
    # ontology records for every resolved compound (all sets) -- retrieval only, no tuning
    (out / "ontology").mkdir(parents=True, exist_ok=True)
    n_ont = 0
    for item in bench.annotatable_items:
        res = pipe.resolve(item.input_id)
        if res.mapping.resolved:
            (out / "ontology" / f"{item.input_id.replace(':', '_')}.json").write_text(json.dumps(pipe.ontology_summary(res.mapping), indent=1, ensure_ascii=False))
            n_ont += 1
    dev_pos = [i for i in bench.primary_positives if assignment[i.input_id] == "dev"]
    dev_neg = [i for i in bench.primary_negatives if assignment[i.input_id] == "dev"]
    rows, meta = build_term_coverage(pipe, dev_pos, dev_neg)
    write_tsv(rows, out / "term_coverage_dev.tsv", COVERAGE_COLUMNS)
    (out / "term_coverage_dev.meta.json").write_text(json.dumps(meta, indent=1))
    print(f"ontology records written: {n_ont}; coverage rows: {len(rows)} -> {out / 'term_coverage_dev.tsv'}")
    print(json.dumps(meta, indent=1))
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
    c = sub.add_parser("split", help="create the deterministic dev/holdout split (output/split.json)")
    c.add_argument("--input", default="data/benchmark/antibiotic_ontology_rule_evaluation.xlsx")
    c.add_argument("--output", default="output")
    c.add_argument("--seed", default="v1")
    c.add_argument("--force", action="store_true")
    c.set_defaults(func=cmd_split)
    d = sub.add_parser("coverage", help="retrieve ontology for resolved compounds and write the DEV term-coverage report")
    d.add_argument("--input", default="data/benchmark/antibiotic_ontology_rule_evaluation.xlsx")
    d.add_argument("--output", default="output")
    d.add_argument("--seed", default="v1")
    d.set_defaults(func=cmd_coverage)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
