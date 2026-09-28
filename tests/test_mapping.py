import pytest

from antibiotic_annotation.cache import JsonFileCache
from antibiotic_annotation.mapping import Mapper, UniChemClient
from antibiotic_annotation.models import ChemicalIdentity
from tests.conftest import FailingTransport


def test_tetracycline_exact_key_with_equivalent_ids(compounds, mapper):
    m = mapper.map(compounds.identity("TAC"))
    assert m.resolved and m.status == "resolved"
    assert m.method == "unichem_inchikey"
    assert m.primary_chebi_id == "CHEBI:27902"
    assert set(m.equivalent_chebi_ids) == {"CHEBI:27902", "CHEBI:77932"}
    assert m.evidence_unioned is True
    assert any("standardised chemical identity" in n for n in m.notes)
    assert m.chebi_ids_for_evidence[0] == "CHEBI:27902"


def test_kirromycin_maps_by_structure_to_mocimycin_entry(compounds, mapper):
    m = mapper.map(compounds.identity("KIR"))
    assert m.method == "unichem_inchikey" and m.primary_chebi_id == "CHEBI:190786"


def test_solithromycin_has_a_chebi_entry(compounds, mapper):
    m = mapper.map(compounds.identity("EM1"))
    assert m.resolved and m.primary_chebi_id == "CHEBI:230261"


def test_component_species_are_preserved(compounds, mapper):
    # neomycin CCD is the single component framycetin (neomycin B), not the mixture
    m = mapper.map(compounds.identity("NMY"))
    assert m.primary_chebi_id == "CHEBI:7508" and "CHEBI:7507" not in m.equivalent_chebi_ids
    # gentamicin CCD is gentamicin C1a, not the gentamicin family
    m = mapper.map(compounds.identity("LLL"))
    assert m.primary_chebi_id == "CHEBI:27784" and "CHEBI:17833" not in m.equivalent_chebi_ids
    # virginiamycin CCD is virginiamycin M1 (pristinamycin IIA), not the mixture
    m = mapper.map(compounds.identity("VIR"))
    assert m.primary_chebi_id == "CHEBI:9997" and "CHEBI:87209" not in m.equivalent_chebi_ids


def test_unverified_crossref_is_not_accepted(compounds, mapper):
    # a record that claims both "water" and "oxygen atom" as ChEBI cross-references: only the
    # entry whose InChIKey equals the deposited one is accepted
    ident = compounds.identity("HOH")
    ident.xrefs["ChEBI"] = ["CHEBI:25805", "CHEBI:15377"]
    m = mapper.map(ident)
    assert m.method == "direct_ccd_crossref" and m.primary_chebi_id == "CHEBI:15377"
    assert "CHEBI:25805" not in m.equivalent_chebi_ids
    rejected = [c for c in m.candidates if c.chebi_id == "CHEBI:25805"]
    assert rejected and rejected[0].accepted is False and rejected[0].mismatches == ["inchikey_not_equal"]
    # PDBe's own record carries a verified ChEBI cross-link for water, so the natural path is direct_ccd_crossref
    m2 = mapper.map(compounds.identity("HOH"))
    assert m2.method == "direct_ccd_crossref" and m2.primary_chebi_id == "CHEBI:15377" and m2.confidence == "high"
    # and a record without any cross-link maps by exact InChIKey
    ident = compounds.identity("HOH")
    ident.xrefs.pop("ChEBI")
    m3 = mapper.map(ident)
    assert m3.method == "unichem_inchikey" and m3.primary_chebi_id == "CHEBI:15377"


def test_protonation_only_difference_is_accepted(compounds, mapper):
    m = mapper.map(compounds.identity("T1C"))  # deposited tigecycline carries +2 charge
    assert m.resolved and m.method == "unichem_connectivity_protonation"
    assert "CHEBI:149836" in m.equivalent_chebi_ids
    for c in m.candidates:
        if c.accepted:
            assert set(c.mismatches) <= {"protonation", "charge", "HAtoms", "isotopicExchangeableH"}


def test_stereo_mismatch_is_an_identity_conflict(compounds, mapper):
    m = mapper.map(compounds.identity("HY0"))  # hygromycin B, stereo layer differs from ChEBI:16976
    assert not m.resolved
    assert m.status == "unresolved_identity_conflict"
    assert m.primary_chebi_id is None
    hyg = [c for c in m.candidates if c.chebi_id == "CHEBI:16976"]
    assert hyg and hyg[0].accepted is False
    assert any(k.lower().startswith("stereo") for k in hyg[0].mismatches)


def test_avilamycin_c_is_never_substituted_by_avilamycin_a(compounds, mapper):
    m = mapper.map(compounds.identity("6UQ"))
    assert m.status == "unresolved_no_chebi"
    assert "CHEBI:85646" not in [c.chebi_id for c in m.candidates]


def test_bird_entity_with_several_stereo_candidates_stays_unresolved(compounds, mapper):
    m = mapper.map(compounds.identity("PRD_000226"))  # viomycin: CHEBI:15782 (+ its 3+ conjugate) vs a 2-star stereoisomer entry
    assert m.status == "unresolved_identity_conflict"
    assert "conflict_reason=stereo_undefined_in_query" in m.notes
    assert any(n.startswith("multiple_stereo_candidates=2") for n in m.notes)
    assert "CHEBI:15782" in [c.chebi_id for c in m.candidates]
    assert m.confidence == "none"


def test_bird_entity_with_single_plausible_candidate_maps_at_medium_confidence(compounds, mapper):
    m = mapper.map(compounds.identity("PRD_000505"))  # quinupristin
    assert m.resolved
    assert m.method == "unichem_connectivity_stereo_undefined"
    assert m.confidence == "medium"
    assert m.primary_chebi_id == "CHEBI:8732"
    m = mapper.map(compounds.identity("PRD_000193"))  # capreomycin IA
    assert m.resolved and m.method == "unichem_connectivity_stereo_undefined" and m.primary_chebi_id == "CHEBI:218527"


def test_stereo_undefined_rule_only_applies_to_stereo_less_queries(compounds, mapper):
    m = mapper.map(compounds.identity("HY0"))  # has a stereo layer that conflicts -> never accepted
    assert m.status == "unresolved_identity_conflict" and m.method == "unresolved"


def test_mapping_confidence_levels(compounds, mapper):
    assert mapper.map(compounds.identity("TAC")).confidence == "high"
    assert mapper.map(compounds.identity("T1C")).confidence == "medium"
    assert mapper.map(compounds.identity("6UQ")).confidence == "none"


def test_chebi_input_bypasses_mapping(compounds, mapper):
    m = mapper.map(compounds.identity("CHEBI:26710"))
    assert m.resolved and m.method == "chebi_input" and m.primary_chebi_id == "CHEBI:26710"


def test_network_failure_is_reported_not_negative(compounds, cache):
    ident = compounds.identity("TAC")
    mapper = Mapper(UniChemClient(FailingTransport(), cache))
    m = mapper.map(ident)
    assert m.status == "unresolved_network" and not m.resolved


def test_manual_synonym_only_without_structural_identifier(cache, transport):
    table = {"XXX": {"chebi_id": "CHEBI:27902", "note": "test entry"}, "TAC": {"chebi_id": "CHEBI:1", "note": "must be ignored"}}
    mapper = Mapper(UniChemClient(transport, cache), synonym_table=table)
    assert mapper.map(ChemicalIdentity(input_id="XXX", entity_kind="ccd_nonpolymer", name="mystery")).confidence == "low"
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


def test_unichem_lookups_are_cached(compounds, transport, cache):
    mapper = Mapper(UniChemClient(transport, cache))
    ident = compounds.identity("KSG")
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
    tac = pipe.resolve("TAC").family.forms
    assert len(tac) == len({(f["chebi_id"], f["relation"]) for f in tac})  # no duplicate form entries


def test_related_parent_is_attached_but_never_used_as_identity(tmp_path, transport):
    import json

    from antibiotic_annotation.pipeline import Pipeline

    cfg = tmp_path / "related.json"
    cfg.write_text(json.dumps({"5I0": {"chebi_id": "CHEBI:17076", "name": "streptomycin", "relation": "hydrated, protonated covalent form of"}}))
    pipe = Pipeline(cache_dir=tmp_path / "cache", transport=transport, synonyms=None, related_parents=cfg)
    res = pipe.resolve("5I0")
    assert res.mapping.status == "unresolved_no_chebi"
    assert res.mapping.primary_chebi_id is None
    assert res.mapping.related_parent["chebi_id"] == "CHEBI:17076"


def test_pipeline_unknown_id(tmp_path, transport):
    from antibiotic_annotation.pipeline import Pipeline

    pipe = Pipeline(cache_dir=tmp_path / "cache", transport=transport, synonyms=None)
    res = pipe.resolve("ZZZZ9")
    assert res.mapping.status == "unresolved_identity_not_found"
