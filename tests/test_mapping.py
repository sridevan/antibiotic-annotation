import pytest

from antibiotic_annotation.cache import JsonFileCache
from antibiotic_annotation.mapping import Mapper, UniChemClient
from antibiotic_annotation.models import ChemicalIdentity
from tests.conftest import FailingTransport


def test_tetracycline_exact_key_with_equivalent_ids(rcsb, mapper):
    m = mapper.map(rcsb.identity("TAC"))
    assert m.resolved and m.status == "resolved"
    assert m.method == "unichem_inchikey"
    assert m.primary_chebi_id == "CHEBI:27902"
    assert set(m.equivalent_chebi_ids) == {"CHEBI:27902", "CHEBI:77932"}
    assert m.evidence_unioned is True
    assert any("standardised chemical identity" in n for n in m.notes)
    assert m.chebi_ids_for_evidence[0] == "CHEBI:27902"


def test_kirromycin_maps_by_structure_to_mocimycin_entry(rcsb, mapper):
    m = mapper.map(rcsb.identity("KIR"))
    assert m.method == "unichem_inchikey" and m.primary_chebi_id == "CHEBI:190786"


def test_solithromycin_has_a_chebi_entry(rcsb, mapper):
    m = mapper.map(rcsb.identity("EM1"))
    assert m.resolved and m.primary_chebi_id == "CHEBI:230261"


def test_component_species_are_preserved(rcsb, mapper):
    # neomycin CCD is the single component framycetin (neomycin B), not the mixture
    m = mapper.map(rcsb.identity("NMY"))
    assert m.primary_chebi_id == "CHEBI:7508" and "CHEBI:7507" not in m.equivalent_chebi_ids
    # gentamicin CCD is gentamicin C1a, not the gentamicin family
    m = mapper.map(rcsb.identity("LLL"))
    assert m.primary_chebi_id == "CHEBI:27784" and "CHEBI:17833" not in m.equivalent_chebi_ids
    # virginiamycin CCD is virginiamycin M1 (pristinamycin IIA), not the mixture
    m = mapper.map(rcsb.identity("VIR"))
    assert m.primary_chebi_id == "CHEBI:9997" and "CHEBI:87209" not in m.equivalent_chebi_ids


def test_unverified_rcsb_crossref_is_not_accepted(rcsb, mapper):
    # RCSB lists both "water" and "oxygen atom" as ChEBI cross-references for HOH
    m = mapper.map(rcsb.identity("HOH"))
    assert m.primary_chebi_id == "CHEBI:15377"
    assert "CHEBI:25805" not in m.equivalent_chebi_ids
    rejected = [c for c in m.candidates if c.chebi_id == "CHEBI:25805"]
    assert rejected and rejected[0].accepted is False


def test_protonation_only_difference_is_accepted(rcsb, mapper):
    m = mapper.map(rcsb.identity("T1C"))  # deposited tigecycline carries +2 charge
    assert m.resolved and m.method == "unichem_connectivity_protonation"
    assert "CHEBI:149836" in m.equivalent_chebi_ids
    for c in m.candidates:
        if c.accepted:
            assert set(c.mismatches) <= {"protonation", "charge", "HAtoms", "isotopicExchangeableH"}


def test_stereo_mismatch_is_an_identity_conflict(rcsb, mapper):
    m = mapper.map(rcsb.identity("HY0"))  # hygromycin B, stereo layer differs from ChEBI:16976
    assert not m.resolved
    assert m.status == "unresolved_identity_conflict"
    assert m.primary_chebi_id is None
    hyg = [c for c in m.candidates if c.chebi_id == "CHEBI:16976"]
    assert hyg and hyg[0].accepted is False
    assert any(k.lower().startswith("stereo") for k in hyg[0].mismatches)


def test_avilamycin_c_is_never_substituted_by_avilamycin_a(rcsb, mapper):
    m = mapper.map(rcsb.identity("6UQ"))
    assert m.status == "unresolved_no_chebi"
    assert "CHEBI:85646" not in [c.chebi_id for c in m.candidates]


def test_bird_entity_without_stereo_layer_is_flagged_not_mapped(rcsb, mapper):
    m = mapper.map(rcsb.identity("PRD_000226"))
    assert m.status == "unresolved_identity_conflict"
    assert "conflict_reason=stereo_undefined_in_query" in m.notes
    assert "CHEBI:15782" in [c.chebi_id for c in m.candidates]


def test_chebi_input_bypasses_mapping(rcsb, mapper):
    m = mapper.map(rcsb.identity("CHEBI:26710"))
    assert m.resolved and m.method == "chebi_input" and m.primary_chebi_id == "CHEBI:26710"


def test_network_failure_is_reported_not_negative(rcsb, cache):
    ident = rcsb.identity("TAC")
    mapper = Mapper(UniChemClient(FailingTransport(), cache))
    m = mapper.map(ident)
    assert m.status == "unresolved_network" and not m.resolved


def test_manual_synonym_only_without_structural_identifier(cache, transport):
    table = {"XXX": {"chebi_id": "CHEBI:27902", "note": "test entry"}, "TAC": {"chebi_id": "CHEBI:1", "note": "must be ignored"}}
    mapper = Mapper(UniChemClient(transport, cache), synonym_table=table)
    no_structure = ChemicalIdentity(input_id="XXX", entity_kind="ccd_nonpolymer", name="mystery")
    m = mapper.map(no_structure)
    assert m.resolved and m.method == "manual_synonym" and m.primary_chebi_id == "CHEBI:27902"
    with_structure = ChemicalIdentity(input_id="TAC", entity_kind="ccd_nonpolymer", inchikey="OFVLGDICTFRJMM-WESIUVDSSA-N")
    m = mapper.map(with_structure)
    assert m.method == "unichem_inchikey" and m.primary_chebi_id == "CHEBI:27902"


def test_no_structure_and_no_synonym_is_unresolved(cache, transport):
    mapper = Mapper(UniChemClient(transport, cache))
    m = mapper.map(ChemicalIdentity(input_id="YYY", entity_kind="ccd_nonpolymer"))
    assert m.status == "unresolved_no_structure"


def test_unichem_lookups_are_cached(rcsb, transport, cache):
    mapper = Mapper(UniChemClient(transport, cache))
    ident = rcsb.identity("KSG")
    mapper.map(ident)
    n = len(transport.calls)
    mapper.map(ident)
    assert len(transport.calls) == n


def test_pipeline_prefers_the_monomeric_specific_entity(tmp_path, transport):
    from antibiotic_annotation.pipeline import Pipeline

    pipe = Pipeline(cache_dir=tmp_path / "cache", transport=transport, synonyms=None)
    res = pipe.resolve("GLC")
    m = res.mapping
    # six ChEBI entries (glucose plus glucan polymer classes) share the InChIKey of the deposited glucose
    assert len(m.equivalent_chebi_ids) >= 3 and m.evidence_unioned
    assert m.primary_chebi_id == "CHEBI:17925"  # alpha-D-glucose, not (1->4)-alpha-D-glucan
    assert res.primary_name == "alpha-D-glucose"


def test_pipeline_reports_family_context_for_components(tmp_path, transport):
    from antibiotic_annotation.pipeline import Pipeline

    pipe = Pipeline(cache_dir=tmp_path / "cache", transport=transport, synonyms=None)
    res = pipe.resolve("NMY")
    assert res.mapping.primary_chebi_id == "CHEBI:7508"
    assert any(f["chebi_id"] == "CHEBI:7507" for f in res.family.member_of)


def test_pipeline_unknown_id(tmp_path, transport):
    from antibiotic_annotation.pipeline import Pipeline

    pipe = Pipeline(cache_dir=tmp_path / "cache", transport=transport, synonyms=None)
    res = pipe.resolve("ZZZZ9")
    assert res.mapping.status == "unresolved_identity_not_found"
