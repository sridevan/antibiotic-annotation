# DEV term-coverage report (summary)

Generated from `term_coverage_dev.tsv` by `scripts/summarise_coverage.py`. ChEBI release 255.
Development split only (`split.json`, seed v1): **12 mapped positives** and **16 mapped negatives** form the denominators.
No holdout or paper compound was inspected. No class roots have been chosen.

Unmapped development positives (not counted in coverage): Evernimycin (EVN: unresolved_identity_conflict); Flopristin (VIF: unresolved_no_chebi); Hygromycin B (HYG: unresolved_identity_conflict); Linopristin (PRD_002101: unresolved_no_chebi); Viomycin (PRD_000226: unresolved_identity_conflict).
Unmapped development negatives: none.

Columns: pos/neg = number of mapped DEV positives/negatives whose evidence reaches the term; depth = minimum / median is_a distance
from the compound (for roles: distance to the term that asserts the role, 0 = asserted on the compound itself);
ChEBI-wide desc. = number of descendants of the term in all of ChEBI (breadth warning); flags = non-antibacterial wording in name or definition.

## 1. Chemical classes with "antibiotic" in the name reached by any DEV compound

| term | name | pos | neg | depth min/med | ChEBI-wide desc. | flags | positive compounds | negative compounds | definition |
|---|---|---|---|---|---|---|---|---|---|
| CHEBI:25105 | macrolide antibiotic | 3 | 0 | 1/1 | 110 | - | Dalfopristin; Erythromycin A; Virginiamycin |  | A macrocyclic lactone with a ring of twelve or more members which exhibits antibiotic activity. |
| CHEBI:23007 | carbohydrate-containing antibiotic | 2 | 0 | 1/1.5 | 111 | - | Clindamycin; Kasugamycin |  | Any carbohydrate derivative that exhibits antibiotic activity. |
| CHEBI:22507 | aminoglycoside antibiotic | 1 | 0 | 1/1 | 86 | - | Kasugamycin |  |  |
| CHEBI:86478 | antibiotic antifungal agent | 1 | 0 | 2/2 | 75 | antifungal | Kasugamycin |  | Heteroorganic entities that are microbial metabolites (or compounds derived from them) which have significant ... |
| CHEBI:87114 | antibiotic fungicide | 1 | 0 | 1/1 | 7 | antifungal | Kasugamycin |  | Any antibiotic antifungal agent that has been used as a fungicide. |

## 2. Family-level chemical classes (no "antibiotic" in the name), >=2 positives, 0 negatives, <=3000 ChEBI-wide descendants

| term | name | pos | neg | depth min/med | ChEBI-wide desc. | flags | positive compounds | negative compounds | definition |
|---|---|---|---|---|---|---|---|---|---|
| CHEBI:26188 | polyketide | 4 | 0 | 2/2.0 | 2552 | - | Dalfopristin; Erythromycin A; Tetracycline; Virginiamycin |  | Natural and synthetic compounds containing alternating carbonyl and methylene groups ('β-polyketones'), biogen... |
| CHEBI:25106 | macrolide | 3 | 0 | 1/1 | 1619 | - | Dalfopristin; Erythromycin A; Virginiamycin |  | A macrocyclic lactone with a ring of twelve or more members derived from a polyketide. |
| CHEBI:35681 | secondary alcohol | 3 | 0 | 1/1 | 2450 | - | Dalfopristin; Spectinomycin; Virginiamycin |  | A secondary alcohol is a compound in which a hydroxy group, ‒OH, is attached to a saturated carbon atom which ... |
| CHEBI:47779 | aminoglycoside | 3 | 0 | 1/1 | 371 | - | CEM-101; Kasugamycin; Neomycin |  |  |
| CHEBI:63944 | macrocyclic lactone | 3 | 0 | 2/2 | 1696 | - | Dalfopristin; Erythromycin A; Virginiamycin |  | Any lactone in which the cyclic carboxylic ester group forms a part of a cyclic macromolecule. |
| CHEBI:140325 | secondary carboxamide | 2 | 0 | 1/1.0 | 2061 | - | Dalfopristin; Virginiamycin |  | A carboxamide resulting from the formal condensation of a carboxylic acid with a primary amine; formula RC(=O)... |
| CHEBI:140326 | tertiary carboxamide | 2 | 0 | 1/1.0 | 255 | - | Dalfopristin; Virginiamycin |  | A carboxamide resulting from the formal condensation of a carboxylic acid with a secondary amine; formula RC(=... |
| CHEBI:35790 | oxazole | 2 | 0 | 2/2.0 | 392 | - | Dalfopristin; Virginiamycin |  | An azole based on a five-membered heterocyclic aromatic skeleton containing one N and one O atom. |
| CHEBI:36683 | organochlorine compound | 2 | 0 | 1/1.0 | 2730 | - | Chloramphenicol; Clindamycin |  | An organochlorine compound is a compound containing at least one carbon-chlorine bond. |
| CHEBI:46812 | 1,3-oxazoles | 2 | 0 | 1/1.0 | 220 | - | Dalfopristin; Virginiamycin |  |  |
| CHEBI:51750 | alpha,beta-unsaturated carboxylic acid amide | 2 | 0 | 2/2.0 | 204 | - | Dalfopristin; Virginiamycin |  | A monocarboxylic amide of general formula R1R2C=CR3‒C(=O)NR4R5  or R1C≡C‒C(=O)NR2R3 in which the amide C=O fun... |
| CHEBI:51751 | enamide | 2 | 0 | 1/1.0 | 198 | - | Dalfopristin; Virginiamycin |  | An α,β-unsaturated carboxylic acid amide of general formula R1R2C=CR3‒C(=O)NR4R5 in which the amide C=O functi... |

## 3. Specific family terms: direct parents (depth 1) covering exactly 1 positive and 0 negatives, <=500 ChEBI-wide descendants

| term | name | pos | neg | depth min/med | ChEBI-wide desc. | flags | positive compounds | negative compounds | definition |
|---|---|---|---|---|---|---|---|---|---|
| CHEBI:139592 | tertiary alpha-hydroxy ketone | 1 | 0 | 1/1 | 230 | - | Tetracycline |  | An α-hydroxy ketone in which the carbonyl group and the hydroxy group are linked by a carbon bearing two organ... |
| CHEBI:144644 | a tetracycline zwitterion | 1 | 0 | 1/1 | 3 | - | Tetracycline |  |  |
| CHEBI:146295 | pyranobenzodioxin | 1 | 0 | 1/1 | 6 | - | Spectinomycin |  | Any organic heterotricyclic compound whose core skeleton consists of a benzodioxin ring that is ortho-fused to... |
| CHEBI:22479 | amino cyclitol glycoside | 1 | 0 | 1/1 | 83 | - | Kasugamycin |  |  |
| CHEBI:23763 | pyrroline | 1 | 0 | 1/1 | 130 | - | Virginiamycin |  | Any organic heteromonocyclic compound with a structure based on a dihydropyrrole. |
| CHEBI:26895 | tetracyclines | 1 | 0 | 1/1 | 98 | - | Tetracycline |  | A subclass of polyketides having an octahydrotetracene-2-carboxamide skeleton, substituted with many hydroxy a... |
| CHEBI:35275 | S-glycosyl compound | 1 | 0 | 1/1 | 103 | - | Clindamycin |  | A glycosyl compound arising formally from the elimination of water from a glycosidic hydroxy group and a S ato... |
| CHEBI:35359 | carboxamidine | 1 | 0 | 1/1 | 73 | - | Kasugamycin |  | Compounds having the structure RC(=NR)NR2. The term is used as a suffix in systematic nomenclature to denote t... |
| CHEBI:35850 | sulfone | 1 | 0 | 1/1 | 145 | - | Dalfopristin |  | An organosulfur compound having the structure RS(=O)2R (R ≠ H). |
| CHEBI:46770 | pyrrolidinecarboxamide | 1 | 0 | 1/1 | 143 | - | Clindamycin |  |  |
| CHEBI:48923 | erythromycin | 1 | 0 | 1/1 | 6 | - | Erythromycin A |  | Any of several wide-spectrum macrolide antibiotics obtained from actinomycete Saccharopolyspora erythraea (for... |
| CHEBI:59770 | cyclic acetal | 1 | 0 | 1/1 | 101 | - | Spectinomycin |  | An acetal in the molecule of which the acetal carbon and one or both oxygen atoms thereon are members of a rin... |
| CHEBI:59780 | cyclic hemiketal | 1 | 0 | 1/1 | 55 | - | Spectinomycin |  | A hemiacetal  having the structure R2C(OH)OR (R ≠ H), derived from a ketone by formal addition of an alcohol t... |
| CHEBI:72588 | semisynthetic derivative | 1 | 0 | 1/1 | 55 | - | Clindamycin |  | Any organic molecular entity derived from a natural product by partial chemical synthesis. |

## 4. Roles covering >=2 positives, or any antibacterial / antimicrobial / antibiotic / antifungal-related role

| term | name | pos | neg | depth min/med | ChEBI-wide desc. | flags | positive compounds | negative compounds | definition |
|---|---|---|---|---|---|---|---|---|---|
| CHEBI:33281 | antimicrobial agent | 9 | 8 | 0/0 | 56 | - | Chloramphenicol; Clindamycin; Dalfopristin; Erythromycin A; Kasugamycin; Neomycin; Spectinomycin; Tetracycline; Virginiamycin | Avibactam; D-glucose; Fluconazole; Formic acid; Glycerol; Hydrogen peroxide; Psilocybin; Xylitol | A substance that kills or slows the growth of microorganisms, including bacteria, viruses, fungi and protozoan... |
| CHEBI:24432 | biological role | 9 | 16 | 0/0 | 1273 | - | Chloramphenicol; Clindamycin; Dalfopristin; Erythromycin A; Kasugamycin; Neomycin; Spectinomycin; Tetracycline; Virginiamycin | ADP; ATP; Avibactam; D-glucose; FAD; Fluconazole; Formic acid; GTP; Glycerol; Heme; Hydrogen peroxide; Ornithine; Phosphate; Psilocybin; Sulfate; Xylitol | A role played by the molecular entity or part thereof within a biological context. |
| CHEBI:50906 | role | 9 | 16 | 0/0 | 1648 | - | Chloramphenicol; Clindamycin; Dalfopristin; Erythromycin A; Kasugamycin; Neomycin; Spectinomycin; Tetracycline; Virginiamycin | ADP; ATP; Avibactam; D-glucose; FAD; Fluconazole; Formic acid; GTP; Glycerol; Heme; Hydrogen peroxide; Ornithine; Phosphate; Psilocybin; Sulfate; Xylitol | A role is particular behaviour which a material entity may exhibit. |
| CHEBI:52206 | biochemical role | 8 | 16 | 0/0.0 | 752 | - | Chloramphenicol; Dalfopristin; Erythromycin A; Kasugamycin; Neomycin; Spectinomycin; Tetracycline; Virginiamycin | ADP; ATP; Avibactam; D-glucose; FAD; Fluconazole; Formic acid; GTP; Glycerol; Heme; Hydrogen peroxide; Ornithine; Phosphate; Psilocybin; Sulfate; Xylitol | A biological role played by the molecular entity or part thereof within a biochemical context. |
| CHEBI:75787 | prokaryotic metabolite | 7 | 6 | 0/0 | 7 | - | Chloramphenicol; Erythromycin A; Kasugamycin; Neomycin; Spectinomycin; Tetracycline; Virginiamycin | D-glucose; FAD; GTP; Glycerol; Hydrogen peroxide; Ornithine | Any metabolite produced during a metabolic reaction in prokaryotes, the taxon that include members of domains ... |
| CHEBI:76969 | bacterial metabolite | 7 | 6 | 0/0 | 5 | - | Chloramphenicol; Erythromycin A; Kasugamycin; Neomycin; Spectinomycin; Tetracycline; Virginiamycin | D-glucose; FAD; GTP; Glycerol; Hydrogen peroxide; Ornithine | Any prokaryotic metabolite produced during a metabolic reaction in bacteria. |
| CHEBI:33232 | application | 7 | 11 | 0/0.0 | 501 | - | Chloramphenicol; Dalfopristin; Kasugamycin; Neomycin; Spectinomycin; Tetracycline; Virginiamycin | ATP; Avibactam; D-glucose; FAD; Fluconazole; Formic acid; Glycerol; Hydrogen peroxide; Ornithine; Psilocybin; Xylitol | Intended use of the molecular entity or part thereof by humans. |
| CHEBI:25212 | metabolite | 7 | 14 | 0/0 | 59 | - | Chloramphenicol; Erythromycin A; Kasugamycin; Neomycin; Spectinomycin; Tetracycline; Virginiamycin | ADP; ATP; D-glucose; FAD; Formic acid; GTP; Glycerol; Heme; Hydrogen peroxide; Ornithine; Phosphate; Psilocybin; Sulfate; Xylitol | Any intermediate or product resulting from metabolism. The term 'metabolite' subsumes the classes commonly kno... |
| CHEBI:36047 | antibacterial drug | 6 | 2 | 0/0.0 | 5 | - | Chloramphenicol; Dalfopristin; Neomycin; Spectinomycin; Tetracycline; Virginiamycin | Avibactam; Psilocybin | A drug used to treat or prevent bacterial infections. |
| CHEBI:36043 | antimicrobial drug | 6 | 3 | 0/0 | 16 | - | Chloramphenicol; Dalfopristin; Neomycin; Spectinomycin; Tetracycline; Virginiamycin | Avibactam; Fluconazole; Psilocybin | A drug used to treat or prevent microbial infections. |
| CHEBI:35441 | antiinfective agent | 6 | 4 | 0/0.0 | 27 | - | Chloramphenicol; Dalfopristin; Neomycin; Spectinomycin; Tetracycline; Virginiamycin | Avibactam; D-glucose; Fluconazole; Psilocybin | A substance used in the prophylaxis or therapy of infectious diseases. |
| CHEBI:33282 | antibacterial agent | 6 | 7 | 0/0 | 7 | - | Chloramphenicol; Dalfopristin; Neomycin; Spectinomycin; Tetracycline; Virginiamycin | Avibactam; D-glucose; Fluconazole; Formic acid; Glycerol; Psilocybin; Xylitol | A substance (or active part thereof) that kills or slows the growth of bacteria. |
| CHEBI:23888 | drug | 6 | 11 | 0/0 | 284 | - | Chloramphenicol; Dalfopristin; Neomycin; Spectinomycin; Tetracycline; Virginiamycin | ATP; Avibactam; D-glucose; FAD; Fluconazole; Formic acid; Glycerol; Hydrogen peroxide; Ornithine; Psilocybin; Xylitol | Any substance which when absorbed into a living organism may modify one or more of its functions. The term is ... |
| CHEBI:52217 | pharmaceutical | 6 | 11 | 0/0 | 295 | - | Chloramphenicol; Dalfopristin; Neomycin; Spectinomycin; Tetracycline; Virginiamycin | ATP; Avibactam; D-glucose; FAD; Fluconazole; Formic acid; Glycerol; Hydrogen peroxide; Ornithine; Psilocybin; Xylitol | Any substance introduced into a living organism with therapeutic or diagnostic purpose. |
| CHEBI:48001 | protein synthesis inhibitor | 4 | 0 | 0/0.0 | 0 | - | Chloramphenicol; Dalfopristin; Kasugamycin; Tetracycline |  | A compound, usually an anti-bacterial agent or a toxin, which inhibits the synthesis of a protein. |
| CHEBI:76932 | pathway inhibitor | 4 | 0 | 0/0.0 | 22 | - | Chloramphenicol; Dalfopristin; Kasugamycin; Tetracycline |  | An enzyme inhibitor that interferes with one or more steps in a metabolic pathway. |
| CHEBI:23924 | enzyme inhibitor | 4 | 2 | 0/0.0 | 604 | - | Chloramphenicol; Dalfopristin; Kasugamycin; Tetracycline | Avibactam; Fluconazole | A compound or agent that combines with an enzyme in such a manner as to prevent the normal substrate-enzyme co... |
| CHEBI:35222 | inhibitor | 4 | 3 | 0/0 | 692 | - | Chloramphenicol; Dalfopristin; Kasugamycin; Tetracycline | Avibactam; Fluconazole; GTP | A substance that diminishes the rate of a chemical reaction. |
| CHEBI:51086 | chemical role | 4 | 7 | 0/0 | 127 | - | Clindamycin; Dalfopristin; Erythromycin A; Spectinomycin | Fluconazole; Formic acid; Glycerol; Hydrogen peroxide; Ornithine; Psilocybin; Xylitol | A role played by the molecular entity or part thereof within a chemical context. |
| CHEBI:66981 | ophthalmology drug | 3 | 1 | 0/0.0 | 1 | - | Chloramphenicol; Neomycin; Tetracycline | Glycerol | Any compound used for the treatment of eye conditions or eye diseases. |
| CHEBI:76971 | Escherichia coli metabolite | 3 | 6 | 0/0 | 0 | - | Chloramphenicol; Neomycin; Tetracycline | D-glucose; FAD; GTP; Glycerol; Hydrogen peroxide; Ornithine | Any bacterial metabolite produced during a metabolic reaction in Escherichia coli. |
| CHEBI:131604 | Mycoplasma genitalium metabolite | 2 | 0 | 0/0.0 | 0 | - | Chloramphenicol; Virginiamycin |  | Any bacterial metabolite produced during a metabolic reaction in Mycoplasma genitalium. |
| CHEBI:35703 | xenobiotic | 2 | 1 | 0/0 | 3 | - | Clindamycin; Erythromycin A | Fluconazole | A xenobiotic (Greek, xenos "foreign"; bios "life") is a compound that is foreign to a living organism. Princip... |
| CHEBI:78298 | environmental contaminant | 2 | 1 | 0/0 | 2 | - | Clindamycin; Erythromycin A | Fluconazole | Any minor or unwanted substance introduced into the environment that can have undesired effects. |
| CHEBI:15339 | acceptor | 2 | 2 | 2/2.0 | 6 | - | Dalfopristin; Spectinomycin | Ornithine; Psilocybin | A molecular entity that can accept an electron, a pair of electrons, an atom or a group from another molecular... |
| CHEBI:22695 | base | 2 | 2 | 2/2.0 | 5 | - | Dalfopristin; Spectinomycin | Ornithine; Psilocybin | A molecular entity having an available pair of electrons capable of forming a covalent bond with a hydron (Brø... |
| CHEBI:35718 | antifungal agent | 2 | 2 | 0/0.5 | 7 | antifungal | Kasugamycin; Tetracycline | D-glucose; Fluconazole | An  antimicrobial agent that destroys fungi by suppressing their ability to grow or reproduce. |
| CHEBI:39142 | Bronsted base | 2 | 2 | 2/2.0 | 2 | - | Dalfopristin; Spectinomycin | Ornithine; Psilocybin | A molecular entity capable of accepting a hydron from a donor (Brønsted acid). |
| CHEBI:86328 | antifungal agrochemical | 1 | 0 | 0/0 | 0 | antifungal | Kasugamycin |  | Any substance used in acriculture, horticulture, forestry, etc. for its fungicidal properties. |
| CHEBI:86327 | antifungal drug | 1 | 1 | 0/0.5 | 0 | antifungal | Tetracycline | Fluconazole | Any antifungal agent used to prevent or treat fungal infections in humans or animals. |

## 5. Generic chemistry classes that score well but are far too broad (>=3 positives, >3000 ChEBI-wide descendants)

Listed so the breadth trap is visible; these are not candidate roots.

| term | name | pos | neg | ChEBI-wide desc. |
|---|---|---|---|---|
| CHEBI:36962 | organochalcogen compound | 12 | 10 | 130110 |
| CHEBI:36963 | organooxygen compound | 12 | 10 | 123765 |
| CHEBI:33285 | heteroorganic entity | 12 | 13 | 159435 |
| CHEBI:33582 | carbon group molecular entity | 12 | 13 | 203518 |
| CHEBI:50860 | organic molecular entity | 12 | 13 | 203378 |
| CHEBI:25806 | oxygen molecular entity | 12 | 14 | 137235 |
| CHEBI:33304 | chalcogen molecular entity | 12 | 14 | 143978 |
| CHEBI:23367 | molecular entity | 12 | 16 | 210611 |
| CHEBI:24431 | chemical entity | 12 | 16 | 216840 |
| CHEBI:33579 | main group molecular entity | 12 | 16 | 208670 |
| CHEBI:33675 | p-block molecular entity | 12 | 16 | 207974 |
| CHEBI:36357 | polyatomic entity | 8 | 16 | 147317 |
| CHEBI:35352 | organonitrogen compound | 7 | 9 | 75434 |
| CHEBI:51143 | nitrogen molecular entity | 7 | 9 | 85407 |
| CHEBI:25367 | molecule | 7 | 10 | 114901 |
| CHEBI:33302 | pnictogen molecular entity | 7 | 10 | 94875 |
| CHEBI:72695 | organic molecule | 7 | 10 | 114349 |
| CHEBI:37577 | heteroatomic molecular entity | 7 | 15 | 63014 |
| CHEBI:33822 | organic hydroxy compound | 6 | 3 | 16936 |
| CHEBI:24651 | hydroxides | 6 | 5 | 28757 |
| CHEBI:33608 | hydrogen molecular entity | 6 | 5 | 31133 |
| CHEBI:33674 | s-block molecular entity | 6 | 5 | 31922 |
| CHEBI:17087 | ketone | 5 | 0 | 9176 |
| CHEBI:63161 | glycosyl compound | 5 | 0 | 9837 |
| CHEBI:36586 | carbonyl compound | 5 | 2 | 41584 |
| CHEBI:36587 | organic oxo compound | 5 | 2 | 41604 |
| CHEBI:63299 | carbohydrate derivative | 5 | 3 | 21762 |
| CHEBI:78616 | carbohydrates and carbohydrate derivatives | 5 | 6 | 29319 |
| CHEBI:24532 | organic heterocyclic compound | 5 | 8 | 61623 |
| CHEBI:33595 | cyclic compound | 5 | 8 | 90337 |
| CHEBI:33832 | organic cyclic compound | 5 | 8 | 89616 |
| CHEBI:5686 | heterocyclic compound | 5 | 8 | 61777 |
| CHEBI:38104 | oxacycle | 4 | 0 | 21836 |
| CHEBI:3992 | cyclic ketone | 4 | 0 | 5110 |
| CHEBI:30879 | alcohol | 4 | 1 | 6046 |
| CHEBI:32988 | amide | 4 | 1 | 40756 |
| CHEBI:33256 | primary amide | 4 | 1 | 38696 |
| CHEBI:37622 | carboxamide | 4 | 1 | 33610 |
| CHEBI:24400 | glycoside | 3 | 0 | 6565 |
| CHEBI:25000 | lactone | 3 | 0 | 6775 |
| CHEBI:33308 | carboxylic ester | 3 | 0 | 16331 |
| CHEBI:51026 | macrocycle | 3 | 0 | 9550 |
| CHEBI:35701 | ester | 3 | 4 | 28603 |
| CHEBI:38101 | organonitrogen heterocyclic compound | 3 | 8 | 38002 |

## 6. Evidence per DEV positive

| compound | id | ChEBI primary | antibiotic-named is_a ancestors | antibacterial-relevant roles |
|---|---|---|---|---|
| Chloramphenicol | CLM | CHEBI:17698 chloramphenicol | none | antibacterial drug; protein synthesis inhibitor; antibacterial agent; antimicrobial agent |
| Clindamycin | CLY | CHEBI:3745 clindamycin | carbohydrate-containing antibiotic (depth 1) | antimicrobial agent |
| Dalfopristin | DOL | CHEBI:4309 dalfopristin | macrolide antibiotic (depth 1) | antibacterial drug; protein synthesis inhibitor; antibacterial agent (via role hierarchy); antimicrobial agent |
| CEM-101 | EM1 | CHEBI:230261 Solithromycin | none | none |
| Erythromycin A | ERY | CHEBI:42355 erythromycin A | macrolide antibiotic (depth 3) | antimicrobial agent |
| Evernimycin | EVN | UNMAPPED (unresolved_identity_conflict) |  |  |
| Hygromycin B | HYG | UNMAPPED (unresolved_identity_conflict) |  |  |
| Kirromycin | KIR | CHEBI:190786 Mocimycin | none | none |
| Kasugamycin | KSG | CHEBI:81419 kasugamycin | aminoglycoside antibiotic (depth 1); antibiotic fungicide (depth 1); carbohydrate-containing antibiotic (depth 2); antibiotic antifungal agent (depth 2) | protein synthesis inhibitor; antimicrobial agent |
| Negamycin | NEG | CHEBI:198310 Negamycin | none | none |
| Neomycin | NMY | CHEBI:7508 framycetin | none | antibacterial drug; antibacterial agent (via role hierarchy); antimicrobial agent (via role hierarchy) |
| Viomycin | PRD_000226 | UNMAPPED (unresolved_identity_conflict) |  |  |
| Linopristin | PRD_002101 | UNMAPPED (unresolved_no_chebi) |  |  |
| Spectinomycin | SCM | CHEBI:9215 spectinomycin | none | antibacterial drug; antibacterial agent (via role hierarchy); antimicrobial agent |
| Tetracycline | TAC | CHEBI:27902 tetracycline | none | antibacterial drug; protein synthesis inhibitor; antibacterial agent; antimicrobial agent |
| Flopristin | VIF | UNMAPPED (unresolved_no_chebi) |  |  |
| Virginiamycin | VIR | CHEBI:9997 pristinamycin IIA | macrolide antibiotic (depth 1) | antibacterial drug; antibacterial agent (via role hierarchy); antimicrobial agent |

## 7. Totals

438 terms in the full report: 283 chemical classes, 155 roles, 0 other.
