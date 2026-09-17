from pathlib import Path

import pytest

from antibiotic_annotation.benchmark import (
    BenchmarkSet,
    load_benchmark,
    stratified_split,
)

WORKBOOK = Path(__file__).resolve().parents[1] / "data" / "benchmark" / "antibiotic_ontology_rule_evaluation.xlsx"


@pytest.fixture(scope="module")
def bench() -> BenchmarkSet:
    return load_benchmark(WORKBOOK)


def names(items):
    return sorted(i.compound for i in items)


def test_primary_positive_set_excludes_peptides_and_ornithine(bench):
    pos = names(bench.primary_positives)
    assert len(pos) == 23
    for excluded in ("Apidaecin", "ErmBL", "SpeFL", "VemP", "Ornithine"):
        assert excluded not in pos
    assert "Tetracycline" in pos and "Viomycin" in pos and "Quinupristin" in pos


def test_prd_entities_keep_entity_kind(bench):
    by_name = {i.compound: i for i in bench.primary_positives}
    assert by_name["Viomycin"].input_id == "PRD_000226"
    assert by_name["Viomycin"].entity_kind == "bird_peptide_like"
    assert by_name["Tetracycline"].input_id == "TAC"
    assert by_name["Tetracycline"].entity_kind == "ccd_nonpolymer"


def test_primary_negatives_are_22_and_include_hard_cases(bench):
    neg = bench.primary_negatives
    assert len(neg) == 22
    ids = {i.input_id for i in neg}
    assert {"HOH", "GTP", "ORN", "TPF", "X8Q", "FYG", "CHEBI:26710"} <= ids
    # fluconazole appears twice in the sheet but is one deposited species
    assert sum(1 for i in neg if i.input_id == "TPF") == 1
    assert all(i.label == 0 for i in neg)
    strata = {i.stratum for i in neg}
    assert {"easy", "hard_ontology", "challenge_negative"} <= strata


def test_clavulanic_acid_is_challenge_not_negative(bench):
    assert names(bench.challenge) == ["Clavulanic acid"]
    assert bench.challenge[0].label is None
    assert "Clavulanic acid" not in names(bench.primary_negatives)


def test_peptide_track_holds_the_four_polymer_peptides(bench):
    assert names(bench.peptide_track) == ["Apidaecin", "ErmBL", "SpeFL", "VemP"]
    by_name = {i.compound: i for i in bench.peptide_track}
    assert by_name["Apidaecin"].label is None  # AMP: separate pathway, no truth here
    assert by_name["ErmBL"].label == 0


def test_paper_external_validation_set(bench):
    paper = bench.paper_positives
    ids = {i.input_id for i in paper}
    # the nine paper-only compounds ...
    assert {"U3B", "YQM", "P8F", "5I0", "AM2", "CA7", "MUL", "G34", "3QB"} <= ids
    # ... plus the two alternative deposited forms of manual positives
    assert {"HY0", "6O1"} <= ids
    # compounds already in the primary set are not duplicated
    assert not ({"TAC", "SCM", "KSG", "LLL", "CLY", "6UQ"} & ids)
    assert len(paper) == 11
    assert all(i.label == 1 and i.benchmark_set == "paper_positive" for i in paper)


def test_workbook_chebi_ids_are_carried(bench):
    by_name = {i.compound: i for i in bench.primary_positives}
    assert by_name["Tetracycline"].workbook_chebi_id == "CHEBI:27902"
    assert by_name["CEM-101"].workbook_chebi_id is None


def test_split_is_deterministic_and_stratified(bench):
    s1 = stratified_split(bench.primary_items, seed="v1", dev_fraction=0.7)
    s2 = stratified_split(bench.primary_items, seed="v1", dev_fraction=0.7)
    assert s1 == s2
    dev = {i.input_id for i in bench.primary_items if s1[i.input_id] == "dev"}
    hold = {i.input_id for i in bench.primary_items if s1[i.input_id] == "holdout"}
    assert dev.isdisjoint(hold) and dev | hold == {i.input_id for i in bench.primary_items}
    n = len(bench.primary_items)
    assert 0.6 * n <= len(dev) <= 0.8 * n
    # holdout contains several antibiotic classes, hard negatives and easy negatives
    hold_items = [i for i in bench.primary_items if i.input_id in hold]
    pos_strata = {i.stratum for i in hold_items if i.label == 1}
    neg_strata = {i.stratum for i in hold_items if i.label == 0}
    assert len(pos_strata) >= 3
    assert {"easy", "hard_ontology"} <= neg_strata
    # a different seed gives a different split
    assert stratified_split(bench.primary_items, seed="v2", dev_fraction=0.7) != s1
