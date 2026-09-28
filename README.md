# antibiotic-annotation

Find the bound chemical entities in a **specific biological assembly** of a PDB entry that have
antibiotic-like properties, and say *why*.

The tool starts from the deposited chemical species (a Chemical Component Dictionary ligand or
a BIRD reference molecule, read from PDBe), maps it to ChEBI by structure through UniChem, walks
the ChEBI ontology, and applies one small, explicit rule. All data comes from EMBL-EBI services
(PDBe, UniChem, OLS/ChEBI). Every decision comes with the ontology
evidence that produced it, and every entity that was inspected but not returned gets a
diagnostic status, so "not an antibiotic" is never confused with "could not be mapped".

Built for ribosome structures; works for any PDB entry.

```
$ antibiotic-annotation 5J7L 1
PDB   Assembly  ID   Type  Name          Antibiotic evidence
5J7L  1         TAC  CCD   tetracycline  antibacterial drug
```

---

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Command line](#command-line)
- [Python API](#python-api)
- [Output format](#output-format)
- [How it works](#how-it-works)
  - [Assembly scope](#1-assembly-scope)
  - [Chemical identity](#2-chemical-identity)
  - [Identity mapping to ChEBI](#3-identity-mapping-to-chebi)
  - [Ontology retrieval](#4-ontology-retrieval)
  - [The antibiotic-like rule](#5-the-antibiotic-like-rule)
  - [Diagnostics](#6-diagnostics)
- [Caching and offline use](#caching-and-offline-use)
- [Known limitations](#known-limitations)
- [Validation benchmark](#validation-benchmark)
- [Development](#development)
- [Project layout](#project-layout)
- [Data sources](#data-sources)

---

## Installation

Requires Python 3.11 or newer.

```bash
git clone git@github.com:sridevan/antibiotic-annotation.git
cd antibiotic-annotation
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[test]"
pytest          # 74 tests, run offline against recorded API responses
```

Runtime dependencies are `requests` and `openpyxl` (the latter only for the validation
benchmark reader).

## Quick start

```python
from antibiotic_annotation import find_antibiotic_entities

for hit in find_antibiotic_entities("5J7L", 1):
    print(hit.entity_id, hit.entity_kind, hit.name, hit.chebi_id, hit.reason)
# TAC CCD tetracycline CHEBI:27902 ['ChEBI antibacterial drug']
```

The first call for a new entry fetches from PDBe, RCSB, UniChem and ChEBI and takes a few
seconds per new compound; every response is cached under `cache/`, so repeat calls are instant.

## Command line

```bash
antibiotic-annotation 5J7L 1                 # antibiotic-like entities in assembly 1 of 5J7L
antibiotic-annotation 5J7L 1 --all           # every inspected entity with its status
antibiotic-annotation 5J7L 1 --json          # full structured result (hits + diagnostics)
antibiotic-annotation --pdb 4V85 --assembly 1
antibiotic-annotation annotate NMY           # one CCD entity
antibiotic-annotation annotate PRD_000226 --kind PRD
antibiotic-annotation --cache /path/to/cache 5J7L 1
python -m antibiotic_annotation.cli 5J7L 1   # same, without the console script
```

`--all` adds a second table listing each inspected entity with its status, ChEBI id and the
naming-stem diagnostic:

```
$ antibiotic-annotation 4V85 1 --all
PDB   Assembly  ID          Type  Name      Antibiotic evidence
4V85  1         PRD_000226  PRD   Viomycin  BIRD:Antibiotic

ID          Type  Name                                      Status           ChEBI        Name stem
GNP         CCD   guanosine 5'-[beta,gamma-imido]triphosph  not_antibiotic   CHEBI:78408  -
MG          CCD   magnesium(2+)                             not_antibiotic   CHEBI:18420  -
PRD_000226  PRD   Viomycin                                  antibiotic_like  -            -mycin
```

Exit code 2 with a message on stderr for an unknown entry or assembly, for example
`assembly '9' not found for 5J7L; available: 1, 2`.

## Python API

```python
from antibiotic_annotation import (
    find_antibiotic_entities,  # (pdb_id, assembly_id) -> list[AntibioticEntity]
    inspect_assembly,          # (pdb_id, assembly_id) -> AssemblyResult (hits + diagnostics)
    annotate_entity,           # (entity_id, entity_kind) -> AnnotationRecord
)
```

| Function | Purpose |
|---|---|
| `find_antibiotic_entities(pdb_id, assembly_id)` | Antibiotic-like CCD and PRD entities present in that assembly. |
| `inspect_assembly(pdb_id, assembly_id)` | Same, plus `entities_inspected` and a `diagnostics` row for every entity. |
| `annotate_entity(entity_id, entity_kind)` | Annotate a single entity outside any structure. `entity_kind` is `"CCD"` or `"PRD"`. |

All three accept `pipeline=Pipeline(cache_dir=...)` to control caching, and every result object
has `.to_dict()` for JSON serialisation.

```python
from antibiotic_annotation import inspect_assembly, annotate_entity

result = inspect_assembly("4V64", 1)
for d in result.diagnostics:
    print(d["entity_id"], d["status"], d["naming_stems"])
# HYG unresolved_identity_conflict ['mycin']
# MG  not_antibiotic []
# ZN  not_antibiotic []

rec = annotate_entity("HYG", "CCD")
for c in rec.mapping["candidates"]:
    print(c["chebi_id"], c["name"], c["accepted"], c["mismatches"])
# CHEBI:16976 hygromycin B False ['stereoSp3']
```

Errors: `EntryNotFound` and `AssemblyNotFound` (from `antibiotic_annotation.assembly`) for bad
ids, `ValueError` for an invalid `entity_kind`, `NetworkUnavailable` when a resource is not
cached and the network cannot be used.

## Output format

Abridged `AntibioticEntity` for tetracycline in 5J7L assembly 1 (`--json` prints the complete
record):

```json
{
  "pdb_id": "5J7L",
  "assembly_id": "1",
  "entity_id": "TAC",
  "entity_kind": "CCD",
  "name": "tetracycline",
  "chebi_id": "CHEBI:27902",
  "antibiotic_like": true,
  "status": "antibiotic_like",
  "reason": ["ChEBI antibacterial drug"],
  "evidence": {
    "antibacterial_drug": true,
    "antibacterial_agent": true,
    "antimicrobial_agent": true,
    "antibiotic_class_ancestors": [],
    "antibiotic_named_ancestors_excluded": [],
    "roles_asserted_on": {"CHEBI:36047": ["CHEBI:27902"], "CHEBI:33281": ["CHEBI:27902", "CHEBI:77932"]},
    "bird_antibiotic": false,
    "bird_class": null
  },
  "mapping": {
    "status": "resolved",
    "method": "unichem_inchikey",
    "confidence": "high",
    "inchikey": "OFVLGDICTFRJMM-WESIUVDSSA-N",
    "equivalent_chebi_ids": ["CHEBI:27902", "CHEBI:77932"],
    "evidence_unioned": true,
    "candidates": [{"chebi_id": "CHEBI:27902", "accepted": true, "mismatches": [], "note": "exact standard InChIKey match"}],
    "related_parent": null,
    "notes": ["2 ChEBI entries share standard InChIKey ...; ontology evidence unioned over all of them"]
  },
  "family": {"member_of": [], "forms": [{"chebi_id": "CHEBI:77932", "name": "tetracycline zwitterion", "relation": "is_tautomer_of"}]},
  "name_flags": {"antibiotic_naming_stem": true, "matches": [{"name": "TETRACYCLINE", "stem": "cycline", "source": "deposited_name", "meaning": "INN stem: tetracyclines"}]},
  "chains": ["CG", "DG"],
  "copies": 2
}
```

Field guide:

| Field | Meaning |
|---|---|
| `entity_id`, `entity_kind` | CCD code (`TAC`) or BIRD id (`PRD_000226`); kind is `CCD` or `PRD`, never conflated. |
| `name` | ChEBI name when mapped, otherwise the deposited name. |
| `chebi_id` | Primary ChEBI id, or `null` when unresolved. |
| `antibiotic_like`, `reason` | The decision and the evidence that caused it. |
| `evidence.antibacterial_drug` | CHEBI:36047 asserted on the entity or inherited from an is_a ancestor. Decisive. |
| `evidence.antibacterial_agent`, `antimicrobial_agent` | CHEBI:33282 / CHEBI:33281. Supporting only, never decisive. |
| `evidence.antibiotic_class_ancestors` | is_a ancestors that are antibiotic chemical classes, with depth. Decisive. |
| `evidence.antibiotic_named_ancestors_excluded` | Antibiotic-named ancestors ignored because their label restricts them to non-antibacterial use (e.g. "antibiotic fungicide"). |
| `evidence.bird_class`, `bird_antibiotic` | BIRD `pdbx_reference_molecule.class`; `true` when it is exactly "Antibiotic". Decisive for PRD entities. |
| `mapping` | How the ChEBI id was obtained: status, method, confidence, every candidate considered with its InChIKey-layer mismatches, and notes. |
| `family.member_of` | Mixtures/families that contain this exact species (e.g. framycetin is part of neomycin). Reported, never substituted. |
| `family.forms` | Conjugate acid/base and tautomer entries of the mapped species. |
| `name_flags` | Diagnostic only: naming stems (-mycin, -micin, -cidin, -cillin, -oxacin, -cycline, -penem, -planin, -bactam) found in names and synonyms. |
| `status` | `antibiotic_like`, `not_antibiotic`, or the mapping status when the entity could not be mapped. |
| `chains`, `copies` | Where the entity sits in the assembly. |

## How it works

```
PDB id + assembly id
  │
  ├─ 1. assembly scope      PDBe assembly + molecules endpoints; PDBe search API for PRD ids
  │                          -> CCD and PRD entities present in that assembly (water excluded)
  │
  └─ for each entity
       ├─ 2. identity         PDBe compound summary -> name, formula, charge, InChI, InChIKey, synonyms,
       │                       cross-links, BIRD class/type for PRD entities
       ├─ 3. mapping          UniChem exact InChIKey, then connectivity search -> ChEBI id(s)
       ├─ 4. ontology         OLS4 (+ ChEBI backend) -> is_a closure, roles, definitions, cached per term
       ├─ 5. rule             antibacterial drug role OR antibiotic chemical class OR BIRD Antibiotic
       └─ 6. diagnostics      status for every entity, naming-stem flag
```

### 1. Assembly scope

Only entities present in the requested biological assembly are inspected, not everything
deposited in the entry.

- PDBe `GET /pdbe/api/pdb/entry/assembly/{pdb_id}` lists each assembly with its entities,
  chains (`in_chains`) and copy numbers.
- PDBe `GET /pdbe/api/pdb/entry/molecules/{pdb_id}` gives the CCD code of each bound entity.
- The PDBe search API (`GET /pdbe/search/pdb/select?q=pdb_id:{pdb_id}`) returns one document per
  entity with `prd_id`, `prd_class`, `prd_name` and `prd_type` for BIRD reference molecules; those
  entity ids are intersected with the assembly's entities.

Example: 4V7T has chloramphenicol in assembly 1 but not in assembly 2; the tool returns it only
for assembly 1.

### 2. Chemical identity

PDBe `GET /pdbe/api/pdb/compound/summary/{id}` serves both CCD codes and PRD ids. For PRD
entries `compound_classes` provides the BIRD class ("antibiotic", "inhibitor", ...) and
`compound_type` the type ("Oligopeptide", ...). Database cross-links (ChEBI, ChEMBL, DrugBank,
PubChem, ...) are kept in the identity record; only a ChEBI cross-link is used, and only after
verification (below).

### 3. Identity mapping to ChEBI

Mapping is structural. Names are never used to map, and a related family member is never
substituted for the deposited species.

| Situation | Status / method | Confidence |
|---|---|---|
| Exact standard InChIKey match in UniChem | `resolved` / `unichem_inchikey` | high |
| The PDBe compound record carries a ChEBI cross-link whose entry has the same InChIKey as the deposited species | `resolved` / `direct_ccd_crossref` | high |
| Only protonation, charge or H-count differ (UniChem connectivity search) | `resolved` / `unichem_connectivity_protonation` | medium |
| PRD descriptor has no stereo layer and exactly one plausible single-species ChEBI identity matches by connectivity | `resolved` / `unichem_connectivity_stereo_undefined` | medium |
| Stereo layer or connectivity differs | `unresolved_identity_conflict` | none |
| Structure known to UniChem but no ChEBI entry | `unresolved_no_chebi` | none |
| No InChIKey in the deposited record | `unresolved_no_structure` | none |
| Network failure | `unresolved_network` | none |

Database cross-links are not identity-verified in general (RCSB's PubChem-assigned ChEBI links,
for example, tie water to "oxygen atom" and glucose to glucan polymer classes), so a ChEBI
cross-link on the compound record is accepted only when that ChEBI entry's InChIKey equals the
deposited one; otherwise it is kept as a rejected candidate and UniChem decides.

Several ChEBI ids can share one standard InChIKey: tautomers (tetracycline and its zwitterion),
protonation forms, and polymer classes whose representative structure is the monomer (glucose
and five glucan classes). They are treated as one *standardised chemical identity for mapping
purposes*: a deterministic primary is chosen (single species before macromolecule/mixture
classes, fewest is_a descendants, highest ChEBI star rating, lowest id) and ontology evidence is
unioned across them; `mapping.evidence_unioned` records that this happened.

The exact deposited species is preserved:

| CCD | Mapped ChEBI entity | Family reported separately |
|---|---|---|
| NMY | framycetin (neomycin B), CHEBI:7508 | neomycin (mixture) |
| LLL | gentamicin C1a, CHEBI:27784 | gentamicin |
| VIR | pristinamycin IIA (virginiamycin M1), CHEBI:9997 | virginiamycin |
| KIR | Mocimycin, CHEBI:190786 (exact structure; ChEBI's name for kirromycin) | |
| T1C | tigecycline(1+) with tigecycline as accepted protonation variant | |

`data/config/related_parents.json` can attach a curated biological parent to a species that is a
distinct covalent form of a known drug (5I0 is hydrated, triply protonated streptomycin) without
claiming chemical identity. `data/config/manual_synonyms.json` is a curated fallback used only
when a record has no structural identifier at all; it is empty.

### 4. Ontology retrieval

Each ChEBI term is fetched once and cached as `cache/chebi/CHEBI_<n>.json`:

- OLS4 v2 (`/ols4/api/v2/ontologies/chebi/classes/{iri}`): label, definition, direct parents,
  precomputed ancestors, `has role` restrictions, synonyms, obsolete flag.
- ChEBI backend (`/chebi/backend/api/public/compound/{n}/`): star rating, InChIKey, and incoming
  structural relations (`has part` from mixtures, conjugate acid/base, tautomer).

The is_a closure is computed locally from cached direct parents (cross-checked against OLS4's
list). Roles are collected from the term and all its is_a ancestors, then closed over the role
hierarchy (antibacterial drug implies antibacterial agent implies antimicrobial agent).

### 5. The antibiotic-like rule

An entity is returned when any of these holds:

1. **Antibacterial drug role.** CHEBI:36047 is asserted on the entity or inherited from an
   is_a ancestor.
2. **Antibiotic chemical class.** An is_a ancestor in the chemical-entity branch of ChEBI whose
   label names an antibiotic class: "macrolide antibiotic", "aminoglycoside antibiotic",
   "carbohydrate-containing antibiotic", "peptide antibiotic", "beta-lactam antibiotic", and so
   on. Labels that restrict the class to non-antibacterial use (antifungal, fungicide,
   antineoplastic, insecticide, ...) are excluded and reported separately.
3. **BIRD class "antibiotic"** (PRD entities only). BIRD's own curation as exposed in PDBe's
   `compound_classes`; recorded as `bird_antibiotic`.

No individual antibiotics are hard-coded. The supporting roles antibacterial agent and
antimicrobial agent are always retrieved and reported but never decide on their own: fluconazole
carries "antibacterial agent" in ChEBI and is correctly not returned.

Why two ChEBI signals: chloramphenicol, spectinomycin and tetracycline have the antibacterial
drug role but sit under generic chemical classes ("tetracyclines", "cyclic ketone"), while
kasugamycin, erythromycin and clindamycin have antibiotic class ancestry but no antibacterial
drug role.

### 6. Diagnostics

`inspect_assembly(...).diagnostics` has one row per inspected entity:

| `status` | Meaning |
|---|---|
| `antibiotic_like` | Returned as a hit. |
| `not_antibiotic` | Mapped to ChEBI; rule not satisfied. |
| `unresolved_identity_conflict` | Structure found but stereo/connectivity disagrees with every ChEBI candidate; candidates retained. |
| `unresolved_no_chebi` | Structure known to UniChem but absent from ChEBI. |
| `unresolved_no_structure`, `unresolved_network`, `unresolved_identity_not_found` | No InChIKey, network failure, unknown id. |

Each row also carries `antibiotic_naming_stem` and `naming_stems`. The stem flag is a review aid,
not evidence: "-mycin" only means "actinomycete product" (mitomycin, rapamycin, natamycin and
nigericin carry it too) and "-bactam" marks beta-lactamase inhibitors, so it never influences
`antibiotic_like`. It exists to surface entities such as kirromycin (`not_antibiotic` because
ChEBI's "Mocimycin" entry has no roles) or hygromycin B (unresolved) for manual review.

## Caching and offline use

- Every HTTP response is stored as JSON under `cache/` (PDBe, UniChem, ChEBI namespaces).
  Override with `--cache DIR` or `Pipeline(cache_dir=DIR)`. The repository ships a populated cache
  for the compounds used during development.
- `ANTIBIOTIC_ANNOTATION_OFFLINE=1` forbids network access; a cache miss raises
  `NetworkUnavailable` instead of fetching.
- The ChEBI release is recorded in `cache/chebi_meta/ontology.json` and on each cached term.
  Delete `cache/chebi/` to refresh against a new ChEBI release.

## Known limitations

- **The rule trusts ChEBI.** CHEBI:36047 currently includes some surprising members (psilocybin,
  sodium chloride) and beta-lactamase inhibitors such as avibactam and clavulanic acid, which are
  antibiotic adjuvants rather than antibiotics. They would be returned as antibiotic-like. No
  guard layer is applied in this version.
- **Known antibiotics without ChEBI evidence** are reported `not_antibiotic`: kirromycin (ChEBI
  "Mocimycin" has no roles), negamycin, solithromycin (CEM-101; 2-star entry with generic
  ancestry). The naming-stem flag helps spot some of these.
- **Stereo disagreements stay unresolved** by design: hygromycin B (HYG, HY0), evernimicin (EVN,
  6O1), retapamulin (G34). Viomycin's PRD descriptor matches two distinct ChEBI stereo entries
  and stays unresolved, but is still returned through its BIRD class.
- **Antibiotics with no ChEBI entry** (avilamycin C, flopristin, linopristin, pentacycline) are
  `unresolved_no_chebi`; PRD entities among them are still returned if BIRD says Antibiotic.
- **Streptomycin 5I0** is a hydrated, triply protonated covalent form whose InChIKey connectivity
  differs from streptomycin; it stays unresolved with a curated `related_parent`.
- **Antimicrobial peptides deposited as plain polymer entities** without a BIRD id (apidaecin,
  and regulatory peptides such as ErmBL, SpeFL, VemP) are not inspected. Polymer entities longer
  than 100 residues are never considered for BIRD lookup.
- Results depend on the live ChEBI release; the cache pins what a given run saw.

## Validation benchmark

`data/benchmark/antibiotic_ontology_rule_evaluation.xlsx` is a manually curated benchmark of
ribosome-bound antibiotics, negative controls and challenge cases. The audited id table
`data/benchmark/benchmark_ids.tsv` maps each compound to the deposited CCD/PRD id (derived from
the entity lists of the PDB entries named in the workbook, not by name matching). It records one
correction: the workbook's CCD for capreomycin, CA7, is the detergent Cymal-7; the deposited
capreomycin is BIRD PRD_000193 (Capreomycin IA).

`antibiotic_annotation.benchmark` reads the workbook and can run the mapping over it;
`tests/test_benchmark.py` checks the set construction (23 primary positives, 22 negatives,
clavulanic acid as a challenge case, peptides on a separate track). The production API never
reads the workbook.

## Development

```bash
pytest -q                               # offline; replays tests/fixtures/http/*.json
python scripts/record_fixtures.py       # re-record fixtures from the live APIs (network)
```

Tests inject a replaying transport (`tests/conftest.py`) fed from recorded real responses,
trimmed to the fields the code uses. Re-record when the APIs or the ChEBI release change, then
review test expectations that depend on ontology content.

## Project layout

```
src/antibiotic_annotation/
  api.py          find_antibiotic_entities / inspect_assembly / annotate_entity, result dataclasses
  assembly.py     PDBe assembly + molecules, PDBe search for PRD ids
  identity.py     PDBe compound summary (CCD and PRD) -> ChemicalIdentity
  mapping.py      UniChem exact / connectivity mapping and acceptance rules
  chebi.py        OLS4 + ChEBI backend client, per-term cache, is_a and role closures
  classifier.py   the antibiotic-like rule, evidence, naming-stem diagnostic
  pipeline.py     wiring, JSON cache, family context
  cli.py          command line
  models.py       shared dataclasses and constants
  cache.py        JSON file cache, HTTP transport with retries and offline switch
  benchmark.py    validation workbook reader
data/config/      related_parents.json, manual_synonyms.json (curated, small)
data/benchmark/   workbook + audited id table
tests/            pytest suite and recorded HTTP fixtures
scripts/          record_fixtures.py
cache/            cached API responses (safe to delete)
```

## Data sources

| Source | Used for |
|---|---|
| [PDBe API](https://www.ebi.ac.uk/pdbe/api/doc/) `pdb/entry/assembly`, `pdb/entry/molecules` | assembly composition, CCD codes of bound entities |
| PDBe API `pdb/compound/summary/{id}` | chemical identity of CCD and PRD entities, BIRD class and type, cross-links |
| [PDBe search API](https://www.ebi.ac.uk/pdbe/api/doc/search.html) `search/pdb/select` | PRD (BIRD) ids of an entry's polymer entities |
| [UniChem](https://www.ebi.ac.uk/unichem/) `api/v1/compounds`, `api/v1/connectivity` | InChIKey to ChEBI mapping with per-layer comparison flags |
| [OLS4](https://www.ebi.ac.uk/ols4/) ChEBI, [ChEBI](https://www.ebi.ac.uk/chebi/) backend API | ontology terms, is_a parents, roles, definitions, star ratings |
| wwPDB [BIRD](https://www.wwpdb.org/data/bird) via PDBe | curated class of peptide-like reference molecules |
