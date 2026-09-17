"""Phase 5: term-coverage report over the DEVELOPMENT split only.

For every ChEBI term reachable from a mapped development compound (is_a ancestors, asserted or
inherited roles, and roles implied by the role hierarchy) report how many development
positives and negatives it covers, so that antibiotic class roots can be chosen from
evidence rather than assumed. Nothing here inspects holdout or paper compounds.
"""
from __future__ import annotations

import re
import statistics
from typing import Any

from .models import BenchmarkItem

ROLE_ROOT = "CHEBI:50906"            # role
CHEMICAL_ENTITY_ROOT = "CHEBI:24431"  # chemical entity

NON_ANTIBACTERIAL_PATTERNS = {
    "antifungal": r"antifungal|fungicid|antimycotic",
    "antineoplastic": r"antineoplastic|antitumou?r|anticancer|antiproliferative",
    "antiviral": r"antiviral|anti-viral|HIV",
    "antiprotozoal": r"antiprotozoal|antimalarial|antiplasmodial|trypanocid|antileishmanial|antitrypanosomal|coccidiostat",
    "antiparasitic": r"anthelmintic|antiparasitic|antinematodal|acaricid|nematicid|insecticid",
    "pesticide_herbicide": r"herbicid|pesticid",
    "immunosuppressive": r"immunosuppress",
}

COVERAGE_COLUMNS = [
    "term_id", "term_name", "term_kind", "antibiotic_in_name", "non_antibacterial_semantics",
    "n_positive", "positive_coverage", "n_negative", "negative_coverage",
    "ancestor_depth_min", "ancestor_depth_median", "reached_via",
    "positive_compounds", "negative_compounds", "stars", "chebi_wide_descendants", "definition",
]


def term_kind(chebi, term_id: str) -> str:
    anc = chebi.ancestors(term_id)
    if term_id == ROLE_ROOT or ROLE_ROOT in anc:
        return "role"
    if term_id == CHEMICAL_ENTITY_ROOT or CHEMICAL_ENTITY_ROOT in anc:
        return "chemical_class"
    return "other"


def semantics_flags(name: str | None, definition: str | None) -> list[str]:
    text = f"{name or ''} {definition or ''}"
    return [k for k, pat in NON_ANTIBACTERIAL_PATTERNS.items() if re.search(pat, text, flags=re.IGNORECASE)]


def compound_evidence(pipeline, mapping) -> dict[str, dict[str, Any]]:
    """{term_id: {"depth": int|None, "via": str}} for one mapped compound.

    Chemical classes: is_a ancestors with their minimal depth. Roles: depth is the is_a distance to
    the term that asserts the role (0 = asserted on the compound itself); roles reached only through
    the role hierarchy get via="role_hierarchy" and the depth of the asserted role they derive from.
    """
    summary = pipeline.ontology_summary(mapping)
    ev: dict[str, dict[str, Any]] = {}
    for a in summary["is_a_ancestors"]:
        ev[a["chebi_id"]] = {"depth": a["min_depth"], "via": "is_a"}
    depth_of_term = {a["chebi_id"]: a["min_depth"] for a in summary["is_a_ancestors"]}
    for cid in summary["chebi_ids"]:
        depth_of_term[cid] = 0
    role_depth: dict[str, int] = {}
    for r, info in summary["inherited_roles"].items():
        d = min(depth_of_term.get(t, 0) for t in info["asserted_on"])
        role_depth[r] = d
        ev[r] = {"depth": d, "via": "has_role" if d == 0 else "has_role_on_ancestor"}
    for r in summary["roles_implied_by_role_hierarchy"]:
        rid = r["chebi_id"]
        # depth of the closest asserted role below it in the role hierarchy
        below = [role_depth[x] for x in role_depth if rid in pipeline.chebi.ancestors(x)]
        ev[rid] = {"depth": min(below) if below else None, "via": "role_hierarchy"}
    return ev


def build_term_coverage(pipeline, positives: list[BenchmarkItem], negatives: list[BenchmarkItem]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Returns (rows, meta). Only mapped compounds count; unmapped are listed in meta."""
    per_compound: dict[str, dict[str, dict[str, Any]]] = {}
    unmapped: dict[str, list[str]] = {"positive": [], "negative": []}
    label_of: dict[str, int] = {}
    for item, lab in [(i, 1) for i in positives] + [(i, 0) for i in negatives]:
        res = pipeline.resolve(item.input_id)
        if not res.mapping.resolved:
            unmapped["positive" if lab else "negative"].append(f"{item.compound} ({item.input_id}: {res.mapping.status})")
            continue
        per_compound[item.compound] = compound_evidence(pipeline, res.mapping)
        label_of[item.compound] = lab
    n_pos = sum(1 for c in label_of.values() if c == 1)
    n_neg = sum(1 for c in label_of.values() if c == 0)
    terms: dict[str, dict[str, Any]] = {}
    for compound, ev in per_compound.items():
        for tid, info in ev.items():
            t = terms.setdefault(tid, {"pos": [], "neg": [], "depths": [], "via": set()})
            (t["pos"] if label_of[compound] else t["neg"]).append(compound)
            if info["depth"] is not None:
                t["depths"].append(info["depth"])
            t["via"].add(info["via"])
    rows: list[dict[str, Any]] = []
    for tid, t in terms.items():
        term = pipeline.chebi.term(tid)
        name = term.name or "?"
        rows.append({
            "term_id": tid,
            "term_name": name,
            "term_kind": term_kind(pipeline.chebi, tid),
            "antibiotic_in_name": "antibiotic" in name.lower(),
            "non_antibacterial_semantics": ",".join(semantics_flags(name, term.definition)),
            "n_positive": len(t["pos"]),
            "positive_coverage": round(len(t["pos"]) / n_pos, 3) if n_pos else 0.0,
            "n_negative": len(t["neg"]),
            "negative_coverage": round(len(t["neg"]) / n_neg, 3) if n_neg else 0.0,
            "ancestor_depth_min": min(t["depths"]) if t["depths"] else "",
            "ancestor_depth_median": statistics.median(t["depths"]) if t["depths"] else "",
            "reached_via": ",".join(sorted(t["via"])),
            "positive_compounds": "; ".join(sorted(t["pos"])),
            "negative_compounds": "; ".join(sorted(t["neg"])),
            "stars": term.stars if term.stars is not None else "",
            "chebi_wide_descendants": term.num_descendants if term.num_descendants is not None else "",
            "definition": (term.definition or "").replace("\n", " "),
        })
    rows.sort(key=lambda r: (-r["n_positive"], r["n_negative"], r["term_kind"], r["term_id"]))
    meta = {"n_dev_positives_mapped": n_pos, "n_dev_negatives_mapped": n_neg, "unmapped": unmapped, "chebi_release": pipeline.chebi.release()}
    return rows, meta
