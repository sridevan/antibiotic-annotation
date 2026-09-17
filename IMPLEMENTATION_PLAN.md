# Implementation plan: PDB CCD -> ChEBI antibacterial-antibiotic annotation

Status: PROPOSAL, awaiting approval. No code written yet.
Date: 2026-09-17. All endpoints below were probed live on this date.

## 0. What the inspection found

### Repository
- No git repository. Folder holds the benchmark workbook, three earlier
  versions of it (`antibiotic_ontology_audit.xlsx`, `..._with_controls.xlsx`,
  `antibiotic_ontology_benchmark.xlsx`; the "(1)" copy is byte-identical),
  the Paternoga et al. 2023 PDF and `Ribosome Structures 2019.xlsx`.
- Python 3.11.9 via pyenv. openpyxl/pandas were not installed (installed
  openpyxl+pandas for inspection only).

### Workbook `antibiotic_ontology_rule_evaluation.xlsx` (sheets used)
- "Manual benchmark audit": 28 rows. 23 small-molecule positives, 1
  antimicrobial peptide (Apidaecin), 3 regulatory peptides (ErmBL, SpeFL,
  VemP), 1 negative (Ornithine). It gives PDB entry IDs but NOT CCD IDs.
- "Negative controls": 13 easy negatives, 8 hard negatives, avibactam
  (challenge negative, truth=0 in "Rule evaluation"), clavulanic acid
  (challenge, excluded from denominator), fluconazole duplicate row.
  Primary precision set = 22 negatives (matches Metrics TN+FP=22).
- "Paper benchmark": 17 compounds with CCD codes; 8 overlap the manual
  set, 9 are new positives (omadacycline U3B, eravacycline YQM,
  pentacycline P8F, streptomycin 5I0, apramycin AM2, capreomycin CA7,
  tiamulin MUL, retapamulin G34, lincomycin 3QB).
- "Rule evaluation"/"Metrics": Rule A/B/C/D preliminary table
  (23 positives, 22 negatives). Will be reproduced programmatically.
- "Ontology terms", "Rule conclusions", "Evaluation plan": constraints
  adopted verbatim (36047 = strong-but-fallible, 33282/33281 = supporting
  only, exclusion layer, 70/30 split, freeze before holdout).

### CCD / PRD identifiers for the manual positives (derived from the
"PDB examples" column through the RCSB entry API, not by name matching)

| Compound (workbook) | ID | Entity kind | Note |
|---|---|---|---|
| Avilamycin C | 6UQ | CCD | no ChEBI entry; avilamycin A (CHEBI:85646) has a DIFFERENT InChIKey connectivity block -> must stay unresolved |
| CEM-101 | EM1 | CCD | maps exactly to CHEBI:230261 "Solithromycin" (2-star). Audit said no ChEBI id; that is now outdated |
| Chloramphenicol | CLM | CCD | CHEBI:17698 |
| Clindamycin | CLY | CCD | CHEBI:3745 |
| Dalfopristin | DOL | CCD | CHEBI:4309 |
| Erythromycin A | ERY | CCD | CHEBI:42355 |
| Evernimycin | EVN (5KCS); paper uses 6O1 | CCD | two CCDs, same connectivity block as CHEBI:29581 but different stereo layers -> stereo mismatch |
| Flopristin | VIF | CCD | no ChEBI entry (PubChem only) |
| Gentamycin | LLL | CCD | exact = CHEBI:27784 gentamycin C1a (component), not CHEBI:17833 mixture |
| Hygromycin B | HYG (4V64); paper uses HY0 | CCD | HYG, HY0 and CHEBI:16976 carry three different stereo layers -> stereo mismatch |
| Kasugamycin | KSG | CCD | CHEBI:81419 |
| Kirromycin | KIR | CCD | exact InChIKey = CHEBI:190786 "Mocimycin"; no synonym step needed |
| Linopristin | PRD_002101 | BIRD peptide-like polymer | not a CCD ligand in 4U1V |
| Negamycin | NEG | CCD | CHEBI:198310 |
| Neomycin | NMY | CCD | exact = CHEBI:7508 framycetin (neomycin B), not CHEBI:7507 mixture |
| Paromomycin | PAR | CCD | CHEBI:7934 (parent, not sulfate) |
| Quinupristin | PRD_000505 | BIRD peptide-like polymer | BIRD class = "Antibiotic" |
| Spectinomycin | SCM | CCD | CHEBI:9215 |
| Telithromycin | TEL | CCD | CHEBI:29688 |
| Tetracycline | TAC | CCD | CHEBI:27902 AND CHEBI:77932 (zwitterion tautomer, same InChIKey) |
| Tigecycline | T1C | CCD | CCD is the +2 protonated form; CHEBI:149836 differs only by protonation layer |
| Viomycin | PRD_000226 | BIRD peptide-like polymer | BIRD class = "Antibiotic"; PRD InChIKey has no stereo layer |
| Virginiamycin | VIR | CCD | exact = CHEBI:9997 pristinamycin IIA (= virginiamycin M1), not CHEBI:87209 mixture |

Apidaecin (5O2R), ErmBL (5JTE), SpeFL (6TBV), VemP are plain polymer
entities with no PRD id: routed to the peptide track / excluded.
Sodium chloride has no single CCD (NA + CL): the control is exercised by
annotating CHEBI:26710 directly, so the CLI must accept ChEBI ids as input.

Paper extras: 5I0 (charge +3, no exact ChEBI), CA7 (CHEBI:3371 has no
structure), P8F and G34 (no ChEBI found by exact key), U3B, YQM, AM2,
MUL, 3QB map exactly.

Multi-hit exact InChIKey cases among negatives: GLC (6 ChEBI ids), GOL (2),
EDO (2), HEM (2). NAD -> CHEBI:44215 NAD zwitterion (workbook says 15846).

## 1. Data sources (verified)

| Purpose | Endpoint | Verified behaviour |
|---|---|---|
| CCD / PRD identity | `GET https://data.rcsb.org/rest/v1/core/chemcomp/{CCD or PRD_xxxxxx}`; batch via `POST https://data.rcsb.org/graphql` `chem_comps(comp_ids:[..])` | name, type, formal charge, formula, InChI, InChIKey, SMILES, synonyms, `rcsb_chem_comp_related` (DrugBank/PubChem/ChEMBL/CCDC; ChEBI never present for our 55 ids), `rcsb_chem_comp_annotation` (WHO ATC codes), and for PRD ids `pdbx_reference_molecule` (class e.g. "Antibiotic", type "Oligopeptide") |
| Entry -> ligand ids (benchmark id table only) | RCSB GraphQL `entries(entry_ids)` nonpolymer_entities / polymer_entities.prd_id | used once to build the id table above |
| InChIKey -> ChEBI | `POST https://www.ebi.ac.uk/unichem/api/v1/compounds` `{"type":"inchikey","compound":KEY}` | returns all sources incl. `chebi` ids; several ChEBI ids possible |
| Relaxed identity (salts/protonation) | `POST https://www.ebi.ac.uk/unichem/api/v1/connectivity` `{"type":"inchikey","compound":KEY,"searchComponents":true,"searchSources":true}` | per-match `comparison` flags (protonation, charge, stereo, HAtoms, connectivity ...) so the mismatch type is explicit |
| UniChem by CCD id | same endpoint, `{"type":"sourceID","compound":"TAC","sourceID":3}` | works but UniChem's PDB record for TAC has a stereo-less key -> used only as cross-check, never as primary. Legacy `unichem/rest/...` is dead (returns HTML) |
| ChEBI term + graph (primary) | `GET https://www.ebi.ac.uk/ols4/api/v2/ontologies/chebi/classes/{double-URL-encoded IRI}` | label, definition, synonyms, `directParent`, `hierarchicalAncestor`, `relatedTo` (RO_0000087 has_role edges, RO_0018036 tautomer etc.), isObsolete, numHierarchicalDescendants. ChEBI release 255, loaded 2026-09-16. ~80 ms/term |
| ChEBI term (secondary / cross-check) | `GET https://www.ebi.ac.uk/chebi/backend/api/public/compound/{num}/`, `.../ontology/parents/{num}/`, `.../ontology/children/{num}/` | stars, default_structure (InChIKey), outgoing (is a, has role, is tautomer of ...) AND incoming relations (mixture/component context), `roles_classification`. Legacy `webservices/chebi/2.0` returns HTTP 500 (dead) |
| Bulk fallback (optional, not default) | `https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi_lite.obo` (55 MB, 2026-09-09) | same `ChebiGraph` interface could be backed by the OBO for whole-PDB runs; deferred |
| PDBe compound API | `https://www.ebi.ac.uk/pdbe/graph-api/compound/summary/{CCD}` | works, no ChEBI cross-links -> not used |

## 2. Data flow

```
input id  (CCD "TAC" | PRD "PRD_000226" | "CHEBI:26710")
  |
  v
identity.py / mapping.py  -------------------------------------------
  RCSB chemcomp -> ChemicalIdentity(ccd_id, name, entity_kind
                    [ccd_nonpolymer | bird_peptide_like], formal_charge,
                    inchi, inchikey, smiles, synonyms, xrefs, atc_codes,
                    bird_class)
  UniChem exact InChIKey -> [ChEBI ids]          method=unichem_inchikey
    none -> UniChem connectivity -> candidates + comparison flags
       only protonation/charge/HAtoms differ -> accept,
                                               method=unichem_connectivity_protonation
       stereo/connectivity differ            -> status=unresolved_stereo_mismatch,
                                               candidates recorded, NOT accepted
    none -> RCSB ChEBI xref if ever present    method=direct_ccd_crossref
    none -> curated synonym table (empty by default, audited)  method=manual_synonym
    none -> status=unresolved_no_chebi
  multi-hit exact key (identical InChIKey => identical entity):
    primary = highest stars, then lowest numeric id; others -> mapping.alternates;
    ontology evidence is the union over exact-key hits.
  |
  v
chebi.py  ------------------------------------------------------------
  ChebiClient (OLS4 primary, ChEBI backend secondary), JSON file cache
    cache/chebi/CHEBI_27902.json  (normalised term: name, definition, stars,
       obsolete, direct parents, direct roles, other relations, incoming
       relations, source payload hashes, retrieved_at, chebi_release)
  ChebiGraph: is_a closure computed locally by recursion over cached
    directParent (cross-checked against OLS4 hierarchicalAncestor);
    direct roles; inherited roles = roles asserted on any is_a ancestor;
    role closure = is_a ancestors of each role term (36047 -> 33282 -> 33281).
    Every referenced id is fetched so every id has a name.
  Offline mode (env ANTIBIOTIC_ANNOTATION_OFFLINE=1): cache only, a miss
    raises OntologyUnavailable -> record status "unresolved_network".
  |
  v
classifier.py  -------------------------------------------------------
  Evidence(class_hits[curated roots matched + full antibiotic-named
          ancestors], roles{antibacterial_drug, antibacterial_agent,
          antimicrobial_agent, mechanism roles}, direct_vs_inherited,
          guards fired, corroboration[ChEBI mechanism role | ATC | BIRD])
  Rules A, B, C, D (+ D2 variant, see 5) -> is_antibiotic, confidence,
  reason string, exclusion_flags.
  |
  v
benchmark.py ---------------------------------------------------------
  workbook reader -> BenchmarkSet(positives, negatives, challenge, peptides)
  id table (data/benchmark/benchmark_ids.tsv, checked in, from section 0)
  deterministic stratified split (dev/holdout), metrics per rule per split,
  mapping coverage vs classifier recall vs end-to-end recall,
  term-coverage report, error analysis TSV, per-compound TSV, JSON records.
```

## 3. Project structure

```
antibiotic_annotation/
  pyproject.toml                 # python>=3.11, requests, openpyxl; pytest, responses (test)
  README.md
  IMPLEMENTATION_PLAN.md
  src/antibiotic_annotation/
    __init__.py
    models.py        # dataclasses: ChemicalIdentity, Mapping, ChebiTerm, Evidence,
                     #   AnnotationRecord (to_json), BenchmarkItem, Metrics
    cache.py         # JsonFileCache(root) + HttpSession (retries, timeouts, offline flag)
    identity.py      # RCSB chemcomp/PRD retrieval  (Phase 3 part 1)
    mapping.py       # UniChem exact/connectivity, disambiguation, salt/protonation policy,
                     #   synonym table, MappingResult with method + status
    chebi.py         # OLS4 + ChEBI backend clients, ChebiGraph (closures), term cache
    rules.py         # curated config loader (roots, guards, corroborators, frozen flag)
    classifier.py    # evidence extraction, rules A-D(+D2), confidence, reasons
    benchmark.py     # workbook reader, id table, split, metrics, reports, error analysis
    coverage.py      # Phase 5 term-coverage report
    cli.py           # annotate / benchmark / coverage / fetch-cache / split
  data/
    benchmark/antibiotic_ontology_rule_evaluation.xlsx   (copied)
    benchmark/benchmark_ids.tsv        # compound -> CCD/PRD id, source PDB entry, note
    config/rules_v1.json               # curated roots + guards (frozen flag), versioned
    config/atc_antibacterial_groups.json  # optional corroborator (see 5)
  cache/
    rcsb/  unichem/  chebi/            # JSON per id, committed for the benchmark
  output/                              # generated (see 8)
  tests/
    conftest.py                        # fake HTTP transport serving tests/fixtures
    fixtures/{rcsb,unichem,chebi}/*.json   # recorded real responses, trimmed
    test_chebi.py  test_mapping.py  test_classifier.py  test_benchmark.py  test_cli.py
  scripts/record_fixtures.py           # one-off: refresh fixtures from live APIs
```

## 4. Phase 5: discovering class roots (data-driven, dev split only)

1. For each DEV positive with a resolved ChEBI id: full is_a ancestor set.
2. Candidate terms = ancestors in the chemical-entity sub-ontology (NOT
   under CHEBI:50906 "role") whose label contains "antibiotic" OR whose
   definition mentions antibiotic/antibacterial, plus family classes that
   the audit already flagged (e.g. CHEBI:26895 tetracyclines) if they
   appear as ancestors.
3. Report per term: id, name, stars, numHierarchicalDescendants,
   n_dev_positives covered + names, n_dev_negatives covered + names,
   definition flags (antifungal / antineoplastic / antiviral wording).
   Preview from OLS4 today (68 labels contain "antibiotic"): heterocyclic
   antibiotic (495 desc), beta-lactam antibiotic (364), carbohydrate-
   containing antibiotic (111), macrolide antibiotic (110), aminoglycoside
   antibiotic (86), peptide antibiotic (67), quinolone (42), sulfonamide
   (24), ... and non-antibacterial branches that must NOT be roots:
   antibiotic antifungal agent/drug (roles), polyene antibiotic
   (antifungals), anthracycline / enediyne antibiotic (antineoplastic).
4. Roots are chosen from that report, written to rules_v1.json with the
   justification, and frozen before the holdout is scored.

## 5. Rules and guards (Phase 6)

- Rule A: 36047 in (direct roles ∪ inherited roles).
- Rule B: any is_a ancestor in curated roots.
- Rule C: A or B.
- Rule D: C and no guard fired. Guards (ontology categories, no compound
  blacklist):
  - G1 inorganic / ionic: is_a descendant of inorganic molecular entity
    (CHEBI:24835) or monoatomic entity -> also `excluded_before_classification`
    together with water (exclusion layer recommended by the workbook).
  - G2 adjuvant: has_role beta-lactamase inhibitor (exact ChEBI id to be
    confirmed from avibactam's record) and no class evidence.
  - G3 role-only-uncorroborated: 36047 present, no class evidence, and no
    corroborating signal. Corroborators, two variants benchmarked:
      D1 (ChEBI only): an antibacterial mechanism role on the compound
         (e.g. protein synthesis inhibitor CHEBI:48001, cell-wall /
         gyrase / RNA-polymerase inhibitor roles; curated list of ChEBI
         role ids in config).
      D2 (ChEBI + WHO ATC via RCSB): D1 OR an ATC code in curated
         antibacterial groups (J01 excluding J01R combinations, A07AA,
         D06A, S01AA ...). Workbook's own "recommended first rule" asks for
         "role + corroborating curated source"; this is the concrete source.
    Expected effect: psilocybin (no mechanism role, no ATC) and avibactam
    fail G3/G2; sodium chloride fails G1. Which positives need D2 rather
    than D1 will be shown as evidence before freezing.
- 33282 / 33281 / BIRD class are recorded as supporting evidence only and
  never make a rule fire.
- Confidence: high = class root AND (36047 or corroborator); medium = one
  signal (class only, or corroborated role only); low = signal present but
  a guard fired or mapping was protonation-relaxed; unknown = unresolved.
  is_antibiotic = Rule D result; the evidence for A/B/C is kept in the record.

## 6. Benchmark, split, metrics (Phases 7-9)

- Positive universe: 23 manual small-molecule positives (primary) + 9
  paper-only positives (secondary). Negatives: 22 unequivocal. Challenge:
  clavulanic acid (reported, not in denominator). Peptides: listed in a
  separate "peptide_track" table with no prediction.
- Split: deterministic (sha256 of the id + fixed seed), stratified by
  (label, stratum) where positive strata = antibiotic class from the paper
  / audit (tetracycline, aminoglycoside, macrolide-ketolide, streptogramin,
  lincosamide, orthosomycin, peptide-like, other) and negative strata =
  easy / hard-ontology / same-domain. ~70/30. Split file written to
  output/split.json and never edited by hand.
- Workflow guard: `benchmark --split holdout` refuses to run unless
  rules_v1.json has `"frozen": true`; coverage report only accepts dev.
- Metrics per rule per split (dev, holdout, all): TP FP TN FN precision
  recall specificity F1, plus mapping coverage, recall among mapped, and
  end-to-end recall. Unresolved mappings are a separate status, never FN.
- Error categories exactly as requested: mapping failure, ontology
  annotation incomplete, class-root missing, role false positive,
  ambiguous biological definition, mixture/salt/component issue; each row
  carries a suggested_action.

## 7. Tests (TDD order)

test_chebi.py: term parse from OLS4 fixture; is_a closure equals fixture
  ancestor list; inherited roles (neomycin-like case where role sits on an
  ancestor); role closure 36047 -> 33282 -> 33281; cache hit makes zero
  HTTP calls; offline miss raises; obsolete term handling; network error
  -> OntologyUnavailable.
test_mapping.py: TAC exact key -> primary 27902, alternate 77932; T1C
  protonation-only difference accepted with method recorded; HY0 stereo
  mismatch -> unresolved with candidates; 6UQ avilamycin C stays
  unresolved (never substituted by avilamycin A); NMY -> framycetin with
  mixture context neomycin; VIR -> pristinamycin IIA; KIR exact -> mocimycin;
  PRD_000226 identity via chemcomp endpoint with entity_kind
  bird_peptide_like; ChEBI-id input bypasses mapping; synonym table only
  used when no structural identifier; network failure -> unresolved_network.
test_classifier.py: TAC (A,B,D true, high); KSG (B true via 22507, A
  false); SCM / spectinomycin (A true, B false, D depends on corroboration
  and the test asserts the reason string); GTP, ORN negative; fluconazole
  33282-only -> all rules false; avibactam A true, D false with adjuvant
  flag; psilocybin A true, D false, reason names G3; sodium chloride
  CHEBI:26710 G1; confidence tiers.
test_benchmark.py: reader excludes Apidaecin/ErmBL/SpeFL/VemP from
  positives, Ornithine is negative, clavulanic acid excluded from
  denominator, 23/22 counts; split deterministic and each stratum present
  in both halves; metrics on synthetic predictions; unresolved not counted
  as FN; error-analysis categorisation.
test_cli.py: `annotate TAC` and `benchmark` run against fixtures/cache in
  offline mode.

## 8. Outputs

output/records/<id>.json (every benchmark compound, target schema from the
brief), output/per_compound.tsv, output/rule_comparison.tsv,
output/metrics.json (dev/holdout/all x rules A,B,C,D1,D2),
output/term_coverage.tsv, output/error_analysis.tsv, output/split.json,
output/peptide_track.tsv.

## 9. Decisions needed before coding

1. Treat viomycin, quinupristin, linopristin as PRD (BIRD) entities inside
   the same pipeline (recommended) rather than deferring them to the
   peptide track.
2. Include the 9 paper-only positives as secondary positives and split
   over the union (recommended).
3. Benchmark both D1 (ChEBI-only corroboration) and D2 (ChEBI + ATC).
4. Relaxed matching policy: accept protonation/charge-only differences,
   reject stereo differences (recommended).
5. `requests` (sync) rather than httpx; openpyxl only, no pandas.
6. `git init` the folder.
