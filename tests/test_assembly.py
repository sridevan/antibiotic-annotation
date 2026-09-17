import pytest

from antibiotic_annotation.assembly import AssemblyClient, AssemblyNotFound, EntryNotFound


def test_assembly_entities_are_scoped_to_the_requested_assembly(transport, cache):
    client = AssemblyClient(transport, cache)
    a1 = {e.entity_id: e for e in client.entities("5J7L", 1)}
    a2 = {e.entity_id: e for e in client.entities("5j7l", "2")}
    assert "TAC" in a1 and "TAC" in a2
    assert a1["TAC"].entity_kind == "CCD" and a1["TAC"].copies == 2 and len(a1["TAC"].chains) == 2
    # additives present only in assembly 1
    assert "MPD" in a1 and "MPD" not in a2
    assert "PUT" in a1 and "PUT" not in a2
    assert "HOH" not in a1
    # entries have two assemblies
    assert client.assembly_ids("5J7L") == ["1", "2"]


def test_prd_entities_are_reported_as_prd(transport, cache):
    ents = {e.entity_id: e for e in AssemblyClient(transport, cache).entities("4V85", 1)}
    assert "PRD_000226" in ents
    assert ents["PRD_000226"].entity_kind == "PRD"
    assert ents["PRD_000226"].name == "Viomycin"
    assert "GNP" in ents and ents["GNP"].entity_kind == "CCD"


def test_invalid_assembly_id_is_a_clear_error(transport, cache):
    with pytest.raises(AssemblyNotFound, match="assembly '9' not found for 5J7L; available: 1, 2"):
        AssemblyClient(transport, cache).entities("5J7L", 9)


def test_unknown_entry_is_a_clear_error(transport, cache):
    with pytest.raises(EntryNotFound):
        AssemblyClient(transport, cache).entities("XXXX", 1)
