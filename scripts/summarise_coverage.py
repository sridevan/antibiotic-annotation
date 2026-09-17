"""Write a Markdown digest of output/term_coverage_dev.tsv (regenerable, no network)."""
from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "output")
rows = list(csv.DictReader((OUT / "term_coverage_dev.tsv").open(encoding="utf-8"), delimiter="\t"))
meta = json.loads((OUT / "term_coverage_dev.meta.json").read_text())
split = json.loads((OUT / "split.json").read_text())
mapping = {r["input_id"]: r for r in csv.DictReader((OUT / "mapping_report.tsv").open(encoding="utf-8"), delimiter="\t")}
npos, nneg = meta["n_dev_positives_mapped"], meta["n_dev_negatives_mapped"]


def table(rs, cols, headers):
    out = ["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)]
    for r in rs:
        out.append("| " + " | ".join(str(c(r)) if callable(c) else str(r[c]) for c in cols) + " |")
    return "\n".join(out) if rs else "_none_"


def defn(r, n=110):
    d = re.sub(r"<[^>]+>", "", r["definition"]).strip()
    return (d[:n] + "...") if len(d) > n else d


common = ["term_id", "term_name", "n_positive", "n_negative", lambda r: f"{r['ancestor_depth_min']}/{r['ancestor_depth_median']}", "chebi_wide_descendants", lambda r: r["non_antibacterial_semantics"] or "-", "positive_compounds", "negative_compounds", defn]
headers = ["term", "name", "pos", "neg", "depth min/med", "ChEBI-wide desc.", "flags", "positive compounds", "negative compounds", "definition"]

cls = [r for r in rows if r["term_kind"] == "chemical_class"]
roles = [r for r in rows if r["term_kind"] == "role"]
A = [r for r in cls if r["antibiotic_in_name"] == "True"]
B = [r for r in cls if r["antibiotic_in_name"] != "True" and int(r["n_positive"]) >= 2 and int(r["n_negative"]) == 0 and int(r["chebi_wide_descendants"] or 0) <= 3000]
C = [r for r in cls if r["antibiotic_in_name"] != "True" and int(r["n_positive"]) == 1 and int(r["n_negative"]) == 0 and r["ancestor_depth_min"] and int(r["ancestor_depth_min"]) <= 1 and int(r["chebi_wide_descendants"] or 0) <= 500]
D = [r for r in roles if int(r["n_positive"]) >= 2 or re.search(r"antibac|antimicrob|antibiot|antiinfect|antitubercul|protein synthesis|antifung", r["term_name"], re.I)]
G = [r for r in cls if int(r["n_positive"]) >= 3 and int(r["chebi_wide_descendants"] or 0) > 3000]

per = []
for iid, part in sorted(split["assignment"].items()):
    it = split["items"][iid]
    if it["label"] != 1 or part != "dev":
        continue
    m = mapping[iid]
    if m["mapping_status"] != "resolved":
        per.append((it["compound"], iid, f"UNMAPPED ({m['mapping_status']})", "", ""))
        continue
    o = json.loads((OUT / "ontology" / f"{iid.replace(':', '_')}.json").read_text())
    named = "; ".join(f"{a['name']} (depth {a['min_depth']})" for a in o["is_a_ancestors"] if "antibiotic" in a["name"].lower()) or "none"
    inh = set(o["inherited_roles"]) | {x["chebi_id"] for x in o["roles_implied_by_role_hierarchy"]}
    r = []
    for cid, nm in [("CHEBI:36047", "antibacterial drug"), ("CHEBI:48001", "protein synthesis inhibitor"), ("CHEBI:33282", "antibacterial agent"), ("CHEBI:33281", "antimicrobial agent")]:
        if cid in inh:
            r.append(nm + ("" if cid in o["inherited_roles"] else " (via role hierarchy)"))
    per.append((it["compound"], iid, m["primary_chebi_id"] + " " + m["primary_name"], named, "; ".join(r) or "none"))

md = f"""# DEV term-coverage report (summary)

Generated from `term_coverage_dev.tsv` by `scripts/summarise_coverage.py`. ChEBI release {meta['chebi_release']}.
Development split only (`split.json`, seed v1): **{npos} mapped positives** and **{nneg} mapped negatives** form the denominators.
No holdout or paper compound was inspected. No class roots have been chosen.

Unmapped development positives (not counted in coverage): {"; ".join(meta['unmapped']['positive']) or 'none'}.
Unmapped development negatives: {"; ".join(meta['unmapped']['negative']) or 'none'}.

Columns: pos/neg = number of mapped DEV positives/negatives whose evidence reaches the term; depth = minimum / median is_a distance
from the compound (for roles: distance to the term that asserts the role, 0 = asserted on the compound itself);
ChEBI-wide desc. = number of descendants of the term in all of ChEBI (breadth warning); flags = non-antibacterial wording in name or definition.

## 1. Chemical classes with "antibiotic" in the name reached by any DEV compound

{table(A, common, headers)}

## 2. Family-level chemical classes (no "antibiotic" in the name), >=2 positives, 0 negatives, <=3000 ChEBI-wide descendants

{table(B, common, headers)}

## 3. Specific family terms: direct parents (depth 1) covering exactly 1 positive and 0 negatives, <=500 ChEBI-wide descendants

{table(C, common, headers)}

## 4. Roles covering >=2 positives, or any antibacterial / antimicrobial / antibiotic / antifungal-related role

{table(D, common, headers)}

## 5. Generic chemistry classes that score well but are far too broad (>=3 positives, >3000 ChEBI-wide descendants)

Listed so the breadth trap is visible; these are not candidate roots.

{table(G, ["term_id", "term_name", "n_positive", "n_negative", "chebi_wide_descendants"], ["term", "name", "pos", "neg", "ChEBI-wide desc."])}

## 6. Evidence per DEV positive

| compound | id | ChEBI primary | antibiotic-named is_a ancestors | antibacterial-relevant roles |
|---|---|---|---|---|
""" + "\n".join(f"| {a} | {b} | {c} | {d} | {e} |" for a, b, c, d, e in per) + f"""

## 7. Totals

{len(rows)} terms in the full report: {len(cls)} chemical classes, {len(roles)} roles, {len(rows) - len(cls) - len(roles)} other.
"""
(OUT / "term_coverage_dev_summary.md").write_text(md, encoding="utf-8")
print(f"wrote {OUT / 'term_coverage_dev_summary.md'} ({len(md.splitlines())} lines)")
