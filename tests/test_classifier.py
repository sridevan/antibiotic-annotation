from antibiotic_annotation.classifier import chebi_evidence, decide, is_antibiotic_class_label, name_stem_flags


def test_antibiotic_class_label_rule():
    assert is_antibiotic_class_label("macrolide antibiotic") == (True, None)
    assert is_antibiotic_class_label("beta-lactam antibiotic") == (True, None)
    assert is_antibiotic_class_label("aliphatic antibiotics") == (True, None)
    ok, why = is_antibiotic_class_label("antibiotic fungicide")
    assert not ok and "antifungal" in why or "fungicid" in why
    assert is_antibiotic_class_label("tetracyclines") == (False, None)
    assert is_antibiotic_class_label(None) == (False, None)


def test_tetracycline_by_role(chebi):
    ev = chebi_evidence(chebi, ["CHEBI:27902", "CHEBI:77932"])
    assert ev.antibacterial_drug and ev.antibacterial_agent and ev.antimicrobial_agent
    assert ev.antibiotic_class_ancestors == []  # ChEBI files it under "tetracyclines", not an antibiotic-named class
    assert "CHEBI:27902" in ev.roles_asserted_on["CHEBI:36047"]
    hit, reasons = decide(ev)
    assert hit and reasons == ["ChEBI antibacterial drug"]


def test_kasugamycin_by_class_ancestry(chebi):
    ev = chebi_evidence(chebi, ["CHEBI:81419"])
    assert not ev.antibacterial_drug
    names = {a["name"] for a in ev.antibiotic_class_ancestors}
    assert "aminoglycoside antibiotic" in names and "carbohydrate-containing antibiotic" in names
    assert {x["name"] for x in ev.antibiotic_named_ancestors_excluded} >= {"antibiotic fungicide"}
    hit, reasons = decide(ev)
    assert hit and reasons[0].startswith("ChEBI antibiotic class: aminoglycoside antibiotic")


def test_spectinomycin_by_role_only(chebi):
    ev = chebi_evidence(chebi, ["CHEBI:9215"])
    assert ev.antibacterial_drug and not ev.antibiotic_class_ancestors
    assert decide(ev)[0]


def test_gtp_is_not_antibiotic(chebi):
    ev = chebi_evidence(chebi, ["CHEBI:15996"])
    assert not decide(ev)[0]
    assert not ev.antibacterial_drug and not ev.antibiotic_class_ancestors


def test_supporting_roles_alone_do_not_qualify(chebi):
    ev = chebi_evidence(chebi, ["CHEBI:46081"])  # fluconazole: antibacterial agent role, antifungal drug
    assert ev.antibacterial_agent and not ev.antibacterial_drug and not ev.antibiotic_class_ancestors
    assert decide(ev) == (False, [])


def test_bird_antibiotic_is_a_reason_on_its_own(chebi):
    from antibiotic_annotation.classifier import Evidence

    ev = Evidence(bird_antibiotic=True, bird_class="Antibiotic")
    assert decide(ev) == (True, ["BIRD class Antibiotic"])


def test_name_stem_flags_are_diagnostic_only():
    f = name_stem_flags({"deposited_name": ["KIRROMYCIN"], "chebi_name": ["Mocimycin"]})
    assert f["antibiotic_naming_stem"] and {m["stem"] for m in f["matches"]} == {"mycin"}
    assert {m["source"] for m in f["matches"]} == {"deposited_name", "chebi_name"}
    assert name_stem_flags({"deposited_name": ["GUANOSINE-5'-TRIPHOSPHATE"]}) == {"antibiotic_naming_stem": False, "matches": []}
    assert name_stem_flags({"deposited_name": ["avibactam"]})["matches"][0]["stem"] == "bactam"
    assert name_stem_flags({"deposited_name": ["tetracycline"]})["matches"][0]["stem"] == "cycline"
    assert name_stem_flags({"deposited_name": ["rapamycin"]})["antibiotic_naming_stem"] is True  # immunosuppressant: why it is only a flag
    assert not name_stem_flags({"deposited_name": ["Mycin"]})["antibiotic_naming_stem"]  # needs a preceding letter
