import json

from antibiotic_annotation import cli
from antibiotic_annotation.pipeline import Pipeline


def _patched_pipeline(monkeypatch, tmp_path, transport):
    def factory(cache_dir="cache"):
        return Pipeline(cache_dir=tmp_path / "cache", transport=transport, synonyms=None, related_parents=None)

    monkeypatch.setattr(cli, "Pipeline", factory)


def test_cli_table(monkeypatch, tmp_path, transport, capsys):
    _patched_pipeline(monkeypatch, tmp_path, transport)
    assert cli.main(["5J7L", "1"]) == 0
    out = capsys.readouterr().out
    assert "TAC" in out and "CCD" in out and "tetracycline" in out and "antibacterial drug" in out
    assert cli.main(["--pdb", "5J7L", "--assembly", "1", "--all"]) == 0
    out = capsys.readouterr().out
    assert "MPD" in out and "not_antibiotic" in out


def test_cli_json_and_errors(monkeypatch, tmp_path, transport, capsys):
    _patched_pipeline(monkeypatch, tmp_path, transport)
    assert cli.main(["4V85", "1", "--json"]) == 0
    d = json.loads(capsys.readouterr().out)
    assert d["hits"][0]["entity_id"] == "PRD_000226" and d["hits"][0]["entity_kind"] == "PRD"
    assert any(x["status"] == "not_antibiotic" for x in d["diagnostics"])
    assert cli.main(["5J7L", "9"]) == 2
    assert "available: 1, 2" in capsys.readouterr().err


def test_cli_annotate(monkeypatch, tmp_path, transport, capsys):
    _patched_pipeline(monkeypatch, tmp_path, transport)
    assert cli.main(["annotate", "TAC"]) == 0
    out = capsys.readouterr().out
    assert "CHEBI:27902" in out and "antibiotic_like: true" in out
    assert cli.main(["annotate", "PRD_000226", "--json"]) == 0
    d = json.loads(capsys.readouterr().out)
    assert d["entity_kind"] == "PRD" and d["evidence"]["bird_antibiotic"] is True
