"""V1 antibiotic-like rule on ChEBI evidence.

An entity is antibiotic-like when its ChEBI entity (or any equivalent id sharing its
standardised identity):

* has the role *antibacterial drug* (CHEBI:36047), asserted on the entity or inherited from an
  is_a ancestor, OR
* has an is_a ancestor that is an antibiotic chemical class: a term of the chemical-entity
  branch whose ChEBI label names it an antibiotic (e.g. "macrolide antibiotic",
  "aminoglycoside antibiotic", "peptide antibiotic", "beta-lactam antibiotic") and whose label
  does not restrict it to non-antibacterial activity (antifungal, antineoplastic, ...).

*antibacterial agent* (CHEBI:33282) and *antimicrobial agent* (CHEBI:33281) are recorded as
supporting evidence only. No individual antibiotics are hard-coded.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

ANTIBACTERIAL_DRUG = "CHEBI:36047"
ANTIBACTERIAL_AGENT = "CHEBI:33282"
ANTIMICROBIAL_AGENT = "CHEBI:33281"
ROLE_ROOT = "CHEBI:50906"

ANTIBIOTIC_LABEL = re.compile(r"\bantibiotics?\b", re.IGNORECASE)
# Labels that name an antibiotic class but explicitly restrict it to a non-antibacterial use.
NON_ANTIBACTERIAL_LABEL = re.compile(r"antifungal|fungicid|antineoplastic|antitumou?r|anticancer|insecticid|acaricid|nematicid|pesticid|herbicid|allergen", re.IGNORECASE)


@dataclass
class Evidence:
    antibacterial_drug: bool = False
    antibacterial_agent: bool = False
    antimicrobial_agent: bool = False
    antibiotic_class_ancestors: list[dict[str, Any]] = field(default_factory=list)
    antibiotic_named_ancestors_excluded: list[dict[str, Any]] = field(default_factory=list)
    roles_asserted_on: dict[str, list[str]] = field(default_factory=dict)  # role id -> ChEBI ids asserting it
    bird_antibiotic: bool = False
    bird_class: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "antibacterial_drug": self.antibacterial_drug,
            "antibacterial_agent": self.antibacterial_agent,
            "antimicrobial_agent": self.antimicrobial_agent,
            "antibiotic_class_ancestors": self.antibiotic_class_ancestors,
            "antibiotic_named_ancestors_excluded": self.antibiotic_named_ancestors_excluded,
            "roles_asserted_on": self.roles_asserted_on,
            "bird_antibiotic": self.bird_antibiotic,
            "bird_class": self.bird_class,
        }


def is_antibiotic_class_label(name: str | None) -> tuple[bool, str | None]:
    """(matches, reason_if_excluded)."""
    if not name or not ANTIBIOTIC_LABEL.search(name):
        return False, None
    m = NON_ANTIBACTERIAL_LABEL.search(name)
    if m:
        return False, f"label restricted to non-antibacterial use ({m.group(0).lower()})"
    return True, None


def chebi_evidence(chebi, chebi_ids: list[str]) -> Evidence:
    """Evidence from ChEBI for one entity; unioned over ids sharing its standardised identity."""
    ev = Evidence()
    seen_anc: dict[str, int] = {}
    for cid in chebi_ids:
        for r, via in chebi.inherited_roles(cid).items():
            ev.roles_asserted_on.setdefault(r, [])
            for v in via:
                if v not in ev.roles_asserted_on[r]:
                    ev.roles_asserted_on[r].append(v)
        for a, d in chebi.ancestor_depths(cid).items():
            if a not in seen_anc or d < seen_anc[a]:
                seen_anc[a] = d
    roles_all = set(ev.roles_asserted_on) | set(chebi.role_closure(ev.roles_asserted_on.keys()))
    ev.antibacterial_drug = ANTIBACTERIAL_DRUG in roles_all
    ev.antibacterial_agent = ANTIBACTERIAL_AGENT in roles_all
    ev.antimicrobial_agent = ANTIMICROBIAL_AGENT in roles_all
    for a, d in sorted(seen_anc.items(), key=lambda kv: (kv[1], kv[0])):
        term = chebi.term(a)
        if not term.name or ROLE_ROOT in chebi.ancestors(a) or a == ROLE_ROOT:
            continue  # roles are handled above; only chemical classes count as class ancestry
        ok, why = is_antibiotic_class_label(term.name)
        row = {"id": a, "name": term.name, "depth": d}
        if ok:
            ev.antibiotic_class_ancestors.append(row)
        elif why:
            ev.antibiotic_named_ancestors_excluded.append({**row, "excluded_because": why})
    return ev


def decide(ev: Evidence) -> tuple[bool, list[str]]:
    reasons: list[str] = []
    if ev.antibacterial_drug:
        reasons.append("ChEBI antibacterial drug")
    if ev.antibiotic_class_ancestors:
        reasons.append("ChEBI antibiotic class: " + "; ".join(f"{a['name']} ({a['id']})" for a in ev.antibiotic_class_ancestors))
    if ev.bird_antibiotic:
        reasons.append("BIRD class Antibiotic")
    return bool(reasons), reasons
