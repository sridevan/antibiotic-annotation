from antibiotic_annotation.coverage import build_term_coverage, semantics_flags, term_kind
from antibiotic_annotation.models import BenchmarkItem
from antibiotic_annotation.pipeline import Pipeline


def _item(name, cid, label):
    return BenchmarkItem(compound=name, label=label, benchmark_set="x", stratum="s", input_id=cid)


def test_semantics_flags():
    assert semantics_flags("antibiotic fungicide", None) == ["antifungal"]
    assert "antineoplastic" in semantics_flags("x", "An antineoplastic agent")
    assert semantics_flags("antibacterial drug", "A drug used to treat bacterial infections.") == []


def test_term_kind(chebi):
    assert term_kind(chebi, "CHEBI:36047") == "role"
    assert term_kind(chebi, "CHEBI:22507") == "chemical_class"


def test_dev_coverage_report(tmp_path, transport):
    pipe = Pipeline(cache_dir=tmp_path / "cache", transport=transport, synonyms=None, related_parents=None)
    pos = [_item("Tetracycline", "TAC", 1), _item("Kasugamycin", "KSG", 1), _item("Spectinomycin", "SCM", 1), _item("Avilamycin C", "6UQ", 1)]
    neg = [_item("GTP", "GTP", 0), _item("Fluconazole", "TPF", 0), _item("Psilocybin", "X8Q", 0)]
    rows, meta = build_term_coverage(pipe, pos, neg)
    assert meta["n_dev_positives_mapped"] == 3 and meta["n_dev_negatives_mapped"] == 3
    assert meta["unmapped"]["positive"] == ["Avilamycin C (6UQ: unresolved_no_chebi)"]
    by_id = {r["term_id"]: r for r in rows}
    agly = by_id["CHEBI:22507"]  # aminoglycoside antibiotic
    assert agly["term_kind"] == "chemical_class" and agly["antibiotic_in_name"] is True
    assert agly["positive_compounds"] == "Kasugamycin" and agly["n_negative"] == 0 and agly["ancestor_depth_min"] == 1
    drug = by_id["CHEBI:36047"]  # antibacterial drug (role)
    assert drug["term_kind"] == "role"
    assert set(drug["positive_compounds"].split("; ")) == {"Spectinomycin", "Tetracycline"}
    assert drug["negative_compounds"] == "Psilocybin"
    agent = by_id["CHEBI:33282"]  # antibacterial agent: reached by fluconazole directly and by others via the role hierarchy
    assert "Fluconazole" in agent["negative_compounds"]
    assert "role_hierarchy" in agent["reached_via"]
    fung = by_id["CHEBI:87114"]  # antibiotic fungicide (kasugamycin is_a)
    assert "antifungal" in fung["non_antibacterial_semantics"]
    assert rows[0]["n_positive"] >= rows[-1]["n_positive"]
