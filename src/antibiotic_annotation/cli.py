"""Command line.

    antibiotic-annotation 5J7L 1              # antibiotic-like entities in assembly 1 of 5J7L
    antibiotic-annotation 5J7L 1 --all        # every inspected entity with its status
    antibiotic-annotation 5J7L 1 --json       # full structured result (hits + diagnostics)
    antibiotic-annotation annotate TAC        # one CCD entity
    antibiotic-annotation annotate PRD_000226 --kind PRD
    antibiotic-annotation --pdb 5J7L --assembly 1
"""
from __future__ import annotations

import argparse
import json
import sys

from .api import annotate_entity, inspect_assembly
from .assembly import AssemblyNotFound, EntryNotFound
from .cache import NetworkUnavailable
from .pipeline import Pipeline


def _evidence_text(rec) -> str:
    parts = []
    ev = rec.evidence
    if ev.antibacterial_drug:
        parts.append("antibacterial drug")
    for a in ev.antibiotic_class_ancestors:
        parts.append(a["name"])
    if ev.bird_antibiotic:
        parts.append("BIRD:Antibiotic")
    if not parts:
        sup = [n for f, n in ((ev.antibacterial_agent, "antibacterial agent"), (ev.antimicrobial_agent, "antimicrobial agent")) if f]
        return ("supporting only: " + ", ".join(sup)) if sup else "-"
    return "; ".join(parts)


def _print_table(rows: list[list[str]], header: list[str]) -> None:
    widths = [max(len(str(r[i])) for r in [header] + rows) for i in range(len(header))]
    fmt = "  ".join("{:<%d}" % w for w in widths)
    print(fmt.format(*header))
    for r in rows:
        print(fmt.format(*[str(x) for x in r]))


def cmd_assembly(args) -> int:
    pipe = Pipeline(cache_dir=args.cache)
    try:
        result = inspect_assembly(args.pdb, args.assembly, pipeline=pipe)
    except (EntryNotFound, AssemblyNotFound, NetworkUnavailable) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps(result.to_dict(), indent=2, ensure_ascii=False))
        return 0
    header = ["PDB", "Assembly", "ID", "Type", "Name", "Antibiotic evidence"]
    rows = [[h.pdb_id, h.assembly_id, h.entity_id, h.entity_kind, (h.name or "")[:40], _evidence_text(h)] for h in result.hits]
    _print_table(rows, header) if rows else print(f"{result.pdb_id} assembly {result.assembly_id}: no antibiotic-like entities among {len(result.entities_inspected)} inspected")
    if args.all:
        print()
        _print_table([[d["entity_id"], d["entity_kind"], (d.get("name") or "")[:40], d["status"], d.get("chebi_id") or "-"] for d in result.diagnostics], ["ID", "Type", "Name", "Status", "ChEBI"])
    return 0


def cmd_annotate(args) -> int:
    pipe = Pipeline(cache_dir=args.cache)
    kind = args.kind or ("PRD" if args.id.upper().startswith("PRD_") else "CCD")
    try:
        rec = annotate_entity(args.id, kind, pipeline=pipe)
    except (ValueError, NetworkUnavailable) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps(rec.to_dict(), indent=2, ensure_ascii=False))
        return 0
    print(f"{rec.entity_id} ({rec.entity_kind})")
    print(rec.name or "(no name)")
    print(rec.chebi_id or "(no ChEBI mapping)")
    print(f"antibiotic_like: {str(rec.antibiotic_like).lower()}")
    print(f"status: {rec.status}")
    print(f"mapping: {rec.mapping['status']} via {rec.mapping['method']} (confidence {rec.mapping['confidence']})")
    print("evidence:")
    ev = rec.evidence
    for a in ev.antibiotic_class_ancestors:
        print(f"  - antibiotic class: {a['name']} ({a['id']}, depth {a['depth']})")
    for flag, label in ((ev.antibacterial_drug, "antibacterial drug"), (ev.antibacterial_agent, "antibacterial agent (supporting)"), (ev.antimicrobial_agent, "antimicrobial agent (supporting)")):
        if flag:
            print(f"  - {label}")
    if ev.bird_class:
        print(f"  - BIRD class: {ev.bird_class}" + (" (antibiotic)" if ev.bird_antibiotic else ""))
    for x in ev.antibiotic_named_ancestors_excluded:
        print(f"  - (excluded) {x['name']}: {x['excluded_because']}")
    for f in rec.family.get("member_of", []):
        print(f"family: component of {f['name']} ({f['chebi_id']})")
    for c in rec.mapping["candidates"]:
        if not c["accepted"]:
            print(f"candidate not accepted: {c['chebi_id']} {c.get('name') or ''} differs in {','.join(c['mismatches'])}")
    for n in rec.mapping["notes"]:
        print(f"note: {n}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="antibiotic-annotation", description="Antibiotic-like CCD/PRD entities in a PDB biological assembly")
    p.add_argument("--cache", default="cache")
    sub = p.add_subparsers(dest="cmd")
    a = sub.add_parser("assembly", help="inspect one biological assembly (default command)")
    a.add_argument("pdb", nargs="?")
    a.add_argument("assembly", nargs="?")
    a.add_argument("--pdb", dest="pdb_flag")
    a.add_argument("--assembly", dest="assembly_flag")
    a.add_argument("--json", action="store_true")
    a.add_argument("--all", action="store_true", help="also list every inspected entity with its status")
    a.set_defaults(func=cmd_assembly)
    b = sub.add_parser("annotate", help="annotate one CCD or PRD entity")
    b.add_argument("id")
    b.add_argument("--kind", choices=["CCD", "PRD"])
    b.add_argument("--json", action="store_true")
    b.set_defaults(func=cmd_annotate)
    return p


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    positional = [x for x in argv if not x.startswith("-")]
    if not positional or positional[0] not in ("assembly", "annotate"):
        argv = ["assembly"] + argv
    args = build_parser().parse_args(argv)
    if args.cmd == "assembly":
        args.pdb = args.pdb or args.pdb_flag
        args.assembly = args.assembly or args.assembly_flag
        if not args.pdb or args.assembly is None:
            build_parser().parse_args(["assembly", "--help"])
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
