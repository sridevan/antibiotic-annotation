# antibiotic-annotation

Find the bound chemical entities in one **biological assembly** of a PDB entry that have
antibiotic-like properties, using the deposited chemical identity (RCSB), UniChem identity
mapping and the ChEBI ontology. Built for ribosome structures; works for any entry.

```python
from antibiotic_annotation import find_antibiotic_entities, inspect_assembly, annotate_entity

hits = find_antibiotic_entities("5J7L", 1)       # -> [AntibioticEntity(entity_id="TAC", entity_kind="CCD", ...)]
result = inspect_assembly("5J7L", 1)             # hits + diagnostics for every inspected entity
record = annotate_entity("PRD_000226", "PRD")    # one CCD or PRD entity
```

```
$ antibiotic-annotation 5J7L 1
PDB   Assembly  ID   Type  Name          Antibiotic evidence
5J7L  1         TAC  CCD   tetracycline  antibacterial drug

$ antibiotic-annotation 4V85 1 --all       # every inspected entity with its status
$ antibiotic-annotation 5J7L 1 --json      # full structured result
$ antibiotic-annotation annotate NMY
$ antibiotic-annotation annotate PRD_000505 --kind PRD
```

## Install

```
pip install -e ".[test]"      # Python >= 3.11; runtime deps: requests, openpyxl (benchmark tests only)
pytest                        # 72 tests, no network (recorded fixtures under tests/fixtures/http)
```

## How it works

```
PDB id + assembly id
  -> PDBe  /pdb/entry/assembly/{id}   entities present in that assembly (chains, copies)
     PDBe  /pdb/entry/molecules/{id}  CCD code of each bound entity
     RCSB  GraphQL polymer_entities    PRD (BIRD) id of short polymer entities  [PDBe exposes no PRD ids]
  -> for each CCD / PRD entity
       RCSB core/chemcomp             deposited name, InChI/InChIKey, BIRD class for PRD entities
       UniChem                        InChIKey -> ChEBI id(s)
       OLS4 + ChEBI backend           is_a parents, roles, definitions (cached per term in cache/chebi/)
  -> V1 rule
       antibiotic_like = has_role antibacterial drug (CHEBI:36047, asserted or inherited)
                      OR is_a ancestor whose ChEBI label names an antibiotic chemical class
                      OR (PRD only) BIRD class == "Antibiotic"
```

Every HTTP response is cached as JSON under `cache/` (override with `--cache` or
`Pipeline(cache_dir=...)`); set `ANTIBIOTIC_ANNOTATION_OFFLINE=1` to forbid network access.

### Identity mapping rules

| Situation | Result |
|---|---|
| exact standard InChIKey match | `resolved`, method `unichem_inchikey` (or `direct_ccd_crossref` when RCSB's PubChem-assigned ChEBI link is verified by InChIKey equality), confidence high |
| only protonation / charge / H-count differ | `resolved`, `unichem_connectivity_protonation`, confidence medium |
| PRD descriptor has no stereo layer and exactly one plausible non-class ChEBI identity matches by connectivity | `resolved`, `unichem_connectivity_stereo_undefined`, confidence medium |
| stereo or connectivity differs | `unresolved_identity_conflict`, candidates and mismatch flags retained |
| no ChEBI entry | `unresolved_no_chebi` |

Several ChEBI ids can share one standard InChIKey (tautomers, zwitterions, polymer classes
built on a monomer). They are treated as one *standardised chemical identity for mapping
purposes*: a deterministic primary is chosen (single species before macromolecule/mixture
classes, fewest descendants, highest star rating, lowest id) and ontology evidence is unioned;
`mapping.evidence_unioned` records this.

The exact deposited species is preserved and never replaced by a family or mixture:
NMY -> framycetin (neomycin B), LLL -> gentamicin C1a, VIR -> pristinamycin IIA (virginiamycin M1).
Mixture membership is reported separately in `family.member_of`. A curated
`data/config/related_parents.json` can attach a biological parent to a distinct covalent form
(5I0 = hydrated, protonated streptomycin) without claiming chemical identity.

### Evidence and diagnostics

`evidence` carries `antibacterial_drug`, the supporting flags `antibacterial_agent` /
`antimicrobial_agent` (never sufficient on their own), `antibiotic_class_ancestors` (with depth),
antibiotic-named ancestors that were excluded because their label restricts them to
non-antibacterial use (e.g. "antibiotic fungicide"), `bird_class` and `bird_antibiotic`.
`reason` lists which evidence caused the match.

`name_flags` is a **diagnostic only**: naming stems such as -mycin, -micin, -cidin, -cillin,
-oxacin, -cycline, -penem, -planin and -bactam found in the deposited name, ChEBI name or their
synonyms, each with its meaning. "-mycin" only means "actinomycete product" (mitomycin,
rapamycin, natamycin and nigericin carry it too) and "-bactam" marks adjuvants, so the flag never
influences `antibiotic_like`; it exists to review `not_antibiotic` and unresolved entities such as
kirromycin or hygromycin B. The `--all` table shows it in the "Name stem" column.

`inspect_assembly(...).diagnostics` has one row per inspected entity with `status`:
`antibiotic_like`, `not_antibiotic` (mapped, rule not met), or the mapping status when the
entity could not be mapped (`unresolved_identity_conflict`, `unresolved_no_chebi`, ...), so
"definitely not" and "could not be mapped" are never confused.

## Known limitations (v1)

* The rule trusts ChEBI. CHEBI:36047 contains a few surprising members (psilocybin, sodium
  chloride) and antibiotic adjuvants such as avibactam and clavulanic acid; these would be
  returned as antibiotic-like. No guard layer is applied in v1.
* Known antibiotics whose ChEBI entry has no antibiotic class and no antibacterial-drug role are
  reported `not_antibiotic` (e.g. kirromycin -> ChEBI "Mocimycin", negamycin, solithromycin).
* Deposited descriptors whose stereo layer disagrees with ChEBI stay unresolved (hygromycin B
  HYG/HY0, evernimicin EVN/6O1, retapamulin G34); viomycin's PRD descriptor matches two ChEBI
  stereo entries and stays unresolved, but is still returned through its BIRD class.
* Antibiotics with no ChEBI entry at all (avilamycin C, flopristin, linopristin, pentacycline)
  are `unresolved_no_chebi`; PRD entities among them are still returned if BIRD says Antibiotic.
* Antimicrobial peptides deposited as plain polymer entities without a BIRD id (e.g. apidaecin)
  are not inspected.

## Benchmark (validation only)

`data/benchmark/antibiotic_ontology_rule_evaluation.xlsx` plus the audited id table
`data/benchmark/benchmark_ids.tsv` (compound -> deposited CCD/PRD id; note the corrected
capreomycin id) are read by `antibiotic_annotation.benchmark` and exercised in
`tests/test_benchmark.py`. The production API never reads the workbook.

## Layout

```
src/antibiotic_annotation/
  api.py          find_antibiotic_entities / inspect_assembly / annotate_entity
  assembly.py     PDBe assembly + molecules, RCSB PRD lookup
  identity.py     RCSB chemcomp (CCD and PRD) -> ChemicalIdentity
  mapping.py      UniChem exact / connectivity mapping rules
  chebi.py        OLS4 + ChEBI backend client, cache, is_a / role closures
  classifier.py   V1 antibiotic-like rule and evidence
  pipeline.py     wiring + JSON cache
  cli.py
  benchmark.py    workbook reader (validation)
scripts/record_fixtures.py   re-record test fixtures from the live APIs
```
