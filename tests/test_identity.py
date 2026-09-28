import pytest

from antibiotic_annotation.cache import NetworkUnavailable
from antibiotic_annotation.identity import IdentityNotFound, CompoundClient, classify_input_id
from tests.conftest import FailingTransport


def test_classify_input_id():
    assert classify_input_id("TAC") == "ccd_nonpolymer"
    assert classify_input_id("PRD_000226") == "bird_peptide_like"
    assert classify_input_id("CHEBI:26710") == "chebi_direct"
    with pytest.raises(ValueError):
        classify_input_id("not an id!")


def test_tetracycline_identity(compounds):
    ident = compounds.identity("TAC")
    assert ident.entity_kind == "ccd_nonpolymer"
    assert ident.name == "TETRACYCLINE"
    assert ident.inchikey == "OFVLGDICTFRJMM-WESIUVDSSA-N"
    assert ident.inchikey_connectivity == "OFVLGDICTFRJMM"
    assert ident.formal_charge == 0
    assert ident.source == "PDBe" and ident.ccd_type == "NON-POLYMER"
    assert "SureChEMBL" in ident.xrefs
    assert any("tetracyclin" in s.lower() for s in ident.synonyms)
    assert ident.bird_class is None and ident.bird_antibiotic is False


def test_bird_entity_keeps_class_and_kind(compounds):
    ident = compounds.identity("PRD_000226")
    assert ident.entity_kind == "bird_peptide_like"
    assert ident.name == "Viomycin"
    assert ident.bird_class == "antibiotic" and ident.bird_antibiotic is True
    assert ident.bird_type == "Oligopeptide"
    assert ident.inchikey.startswith("GXFAIFRPOKBQRV-")


def test_chebi_input_needs_no_network(transport, compounds):
    ident = compounds.identity("CHEBI:26710")
    assert ident.entity_kind == "chebi_direct"
    assert ident.input_id == "CHEBI:26710"
    assert transport.calls == []


def test_unknown_id_raises(compounds):
    with pytest.raises(IdentityNotFound):
        compounds.identity("ZZZZ9")


def test_identity_is_cached(transport, cache):
    client = CompoundClient(transport, cache)
    client.identity("TAC")
    n = len(transport.calls)
    client.identity("TAC")
    assert len(transport.calls) == n
    # a second client with the same cache and a dead network still works
    offline = CompoundClient(FailingTransport(), cache)
    assert offline.identity("TAC").inchikey == "OFVLGDICTFRJMM-WESIUVDSSA-N"


def test_network_failure_propagates(cache):
    client = CompoundClient(FailingTransport(), cache)
    with pytest.raises(NetworkUnavailable):
        client.identity("TAC")
