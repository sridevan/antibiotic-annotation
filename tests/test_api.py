import pytest

from antibiotic_annotation.api import annotate_entity, find_antibiotic_entities, inspect_assembly
from antibiotic_annotation.assembly import AssemblyNotFound, EntryNotFound
from antibiotic_annotation.pipeline import Pipeline


@pytest.fixture
def pipe(tmp_path, transport):
    return Pipeline(cache_dir=tmp_path / "cache", transport=transport, synonyms=None, related_parents=None)


def ids(hits):
    return [h.entity_id for h in hits]


def test_5j7l_assembly_1_returns_tetracycline(pipe):
    hits = find_antibiotic_entities("5J7L", 1, pipeline=pipe)
    assert ids(hits) == ["TAC"]
    tac = hits[0]
    assert tac.pdb_id == "5J7L" and tac.assembly_id == "1"
    assert tac.entity_kind == "CCD" and tac.chebi_id == "CHEBI:27902" and tac.name == "tetracycline"
    assert tac.antibiotic_like is True
    assert tac.evidence.antibacterial_drug and tac.evidence.antibacterial_agent and tac.evidence.antimicrobial_agent
    assert tac.evidence.antibiotic_class_ancestors == [] and tac.evidence.bird_antibiotic is False
    assert tac.reason == ["ChEBI antibacterial drug"]
    assert tac.mapping["status"] == "resolved" and tac.mapping["method"] == "unichem_inchikey"
    assert tac.copies == 2
    d = tac.to_dict()
    assert d["pdb_id"] == "5J7L" and d["evidence"]["antibacterial_drug"] is True


def test_chloramphenicol_and_erythromycin(pipe):
    assert ids(find_antibiotic_entities("4V7T", 1, pipeline=pipe)) == ["CLM"]
    ery = find_antibiotic_entities("4V7U", 1, pipeline=pipe)
    assert ids(ery) == ["ERY"]
    assert any(a["name"] == "macrolide antibiotic" for a in ery[0].evidence.antibiotic_class_ancestors)


def test_ordinary_nucleotides_are_not_returned(pipe):
    r = inspect_assembly("5AFI", 1, pipeline=pipe)  # GDP, fMet, kirromycin, ions
    assert "GDP" in [e.entity_id for e in r.entities_inspected]
    assert "GDP" not in ids(r.hits)
    assert next(d for d in r.diagnostics if d["entity_id"] == "GDP")["status"] == "not_antibiotic"
    r = inspect_assembly("1ATP", 1, pipeline=pipe)
    assert "ATP" in [e.entity_id for e in r.entities_inspected] and "ATP" not in ids(r.hits)
    r = inspect_assembly("4V85", 1, pipeline=pipe)
    assert "GNP" not in ids(r.hits)


def test_prd_antibiotic_returned_with_entity_kind_prd(pipe):
    r = inspect_assembly("4V85", 1, pipeline=pipe)  # viomycin: ChEBI mapping unresolved, BIRD class Antibiotic
    assert ids(r.hits) == ["PRD_000226"]
    v = r.hits[0]
    assert v.entity_kind == "PRD" and v.evidence.bird_antibiotic and v.evidence.bird_class == "Antibiotic"
    assert v.reason == ["BIRD class Antibiotic"]
    assert v.mapping["status"] == "unresolved_identity_conflict" and v.chebi_id is None
    r = inspect_assembly("4U1U", 1, pipeline=pipe)  # quinupristin: mapped (stereo-undefined) + BIRD
    q = next(h for h in r.hits if h.entity_id == "PRD_000505")
    assert q.entity_kind == "PRD" and q.chebi_id == "CHEBI:8732" and "BIRD class Antibiotic" in q.reason
    assert q.mapping["method"] == "unichem_connectivity_stereo_undefined" and q.mapping["confidence"] == "medium"


def test_assembly_specificity(pipe):
    # 4V7T: chloramphenicol is bound only in assembly 1
    assert ids(find_antibiotic_entities("4V7T", 1, pipeline=pipe)) == ["CLM"]
    r2 = inspect_assembly("4V7T", 2, pipeline=pipe)
    assert ids(r2.hits) == [] and "CLM" not in [e.entity_id for e in r2.entities_inspected]
    # 5J7L: MPD is present in assembly 1 only; assembly 2 still contains TAC
    r1, r2 = inspect_assembly("5J7L", 1, pipeline=pipe), inspect_assembly("5J7L", 2, pipeline=pipe)
    assert "MPD" in [e.entity_id for e in r1.entities_inspected] and "MPD" not in [e.entity_id for e in r2.entities_inspected]
    assert ids(r2.hits) == ["TAC"]


def test_invalid_assembly_and_entry(pipe):
    with pytest.raises(AssemblyNotFound, match="available: 1, 2"):
        find_antibiotic_entities("5J7L", 9, pipeline=pipe)
    with pytest.raises(EntryNotFound):
        find_antibiotic_entities("XXXX", 1, pipeline=pipe)


def test_unresolved_mapping_is_a_diagnostic_not_discarded(pipe):
    r = inspect_assembly("4V64", 1, pipeline=pipe)  # hygromycin B: stereo layer conflicts with ChEBI
    assert "HYG" in [e.entity_id for e in r.entities_inspected]
    assert "HYG" not in ids(r.hits)
    d = next(d for d in r.diagnostics if d["entity_id"] == "HYG")
    assert d["status"] == "unresolved_identity_conflict" and d["chebi_id"] is None
    rec = annotate_entity("HYG", "CCD", pipeline=pipe)
    assert rec.status == "unresolved_identity_conflict" and not rec.antibiotic_like
    assert any(c["chebi_id"] == "CHEBI:16976" and not c["accepted"] for c in rec.mapping["candidates"])


def test_known_antibiotic_without_chebi_evidence_is_reported_not_antibiotic(pipe):
    # kirromycin maps exactly to ChEBI "Mocimycin", which carries no roles and no antibiotic class: a documented limitation
    r = inspect_assembly("5AFI", 1, pipeline=pipe)
    d = next(d for d in r.diagnostics if d["entity_id"] == "KIR")
    assert d["status"] == "not_antibiotic" and d["chebi_id"] == "CHEBI:190786"


def test_annotate_entity_validates_kind(pipe):
    with pytest.raises(ValueError):
        annotate_entity("TAC", "XYZ", pipeline=pipe)
    with pytest.raises(ValueError):
        annotate_entity("PRD_000226", "CCD", pipeline=pipe)
    with pytest.raises(ValueError):
        annotate_entity("TAC", "PRD", pipeline=pipe)
    rec = annotate_entity("nmy", "ccd", pipeline=pipe)
    assert rec.entity_id == "NMY" and rec.chebi_id == "CHEBI:7508" and rec.name == "framycetin"
    assert any(f["chebi_id"] == "CHEBI:7507" for f in rec.family["member_of"])
