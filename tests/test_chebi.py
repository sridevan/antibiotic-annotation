import pytest

from antibiotic_annotation.cache import NetworkUnavailable
from antibiotic_annotation.chebi import ChebiClient, ols4_term_url
from tests.conftest import FailingTransport


def test_ols4_url_double_encodes_iri():
    assert ols4_term_url("CHEBI:27902").endswith("http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FCHEBI_27902")


def test_tetracycline_term_merges_both_sources(chebi):
    t = chebi.term("CHEBI:27902")
    assert t.name == "tetracycline"
    assert "antibiotic" in t.definition.lower()
    assert t.stars == 3
    assert set(t.is_a) == {"CHEBI:139592", "CHEBI:26895"}
    assert {"CHEBI:36047", "CHEBI:33282", "CHEBI:33281", "CHEBI:48001"} <= set(t.roles)
    assert ("is_tautomer_of", "CHEBI:77932") in {(r.relation, r.target_id) for r in t.relations}
    assert t.inchikey == "OFVLGDICTFRJMM-WESIUVDSSA-N"
    assert t.chebi_release not in (None, "unknown")


def test_is_a_closure_matches_ols4_precomputed_ancestors(chebi):
    anc = chebi.ancestors("CHEBI:27902")
    assert anc == set(chebi.term("CHEBI:27902").ols4_ancestors)
    assert {"CHEBI:26895", "CHEBI:26188", "CHEBI:24431"} <= anc  # tetracyclines, polyketide, chemical entity
    assert "CHEBI:27902" not in anc
    depths = chebi.ancestor_depths("CHEBI:27902")
    assert depths["CHEBI:26895"] == 1 and depths["CHEBI:24431"] > 1


def test_kasugamycin_class_ancestry_without_drug_role(chebi):
    assert "CHEBI:22507" in chebi.ancestors("CHEBI:81419")  # aminoglycoside antibiotic
    assert chebi.is_a_path("CHEBI:81419", "CHEBI:22507") == ["CHEBI:81419", "CHEBI:22507"]
    assert "CHEBI:48001" in chebi.direct_roles("CHEBI:81419")
    assert "CHEBI:36047" not in chebi.inherited_roles("CHEBI:81419")


def test_spectinomycin_role_without_antibiotic_class(chebi):
    assert "CHEBI:36047" in chebi.direct_roles("CHEBI:9215")
    assert "CHEBI:22507" not in chebi.ancestors("CHEBI:9215")


def test_role_closure_generalises_antibacterial_drug(chebi):
    closure = chebi.role_closure(["CHEBI:36047"])
    assert {"CHEBI:36047", "CHEBI:33282", "CHEBI:33281"} <= closure


def test_inherited_roles_record_the_asserting_term(chebi):
    roles = chebi.inherited_roles("CHEBI:7508")  # framycetin
    assert roles, "framycetin should have roles"
    assert all("CHEBI:7508" in via or via for via in roles.values())
    # every asserting term is the compound itself or one of its ancestors
    anc = chebi.ancestors("CHEBI:7508") | {"CHEBI:7508"}
    assert all(set(v) <= anc for v in roles.values())


def test_mixture_context_via_incoming_relations(chebi):
    inc = {(r.relation, r.target_id) for r in chebi.term("CHEBI:7508").incoming}
    assert ("has_part", "CHEBI:7507") in inc  # neomycin (mixture) has_part framycetin


def test_missing_term_is_handled(chebi):
    assert chebi.has_term("CHEBI:9999999") is False
    assert chebi.ancestors("CHEBI:9999999") == frozenset()


def test_cache_makes_second_client_offline_capable(transport, cache):
    c1 = ChebiClient(transport, cache)
    c1.ancestors("CHEBI:81419")
    c1.role_closure(c1.direct_roles("CHEBI:81419"))
    c2 = ChebiClient(FailingTransport(), cache)
    assert "CHEBI:22507" in c2.ancestors("CHEBI:81419")
    assert c2.term("CHEBI:22507").name == "aminoglycoside antibiotic"
    assert cache.path("chebi", "CHEBI:81419").name == "CHEBI_81419.json"


def test_network_failure_raises(cache):
    c = ChebiClient(FailingTransport(), cache)
    with pytest.raises(NetworkUnavailable):
        c.term("CHEBI:27902")
