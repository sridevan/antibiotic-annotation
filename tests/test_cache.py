import json

import pytest

from antibiotic_annotation.cache import OFFLINE_ENV, JsonFileCache, NetworkUnavailable


def test_cache_roundtrip_and_key_sanitising(tmp_path):
    cache = JsonFileCache(tmp_path)
    cache.put("chebi", "CHEBI:27902", {"name": "tetracycline"})
    assert cache.path("chebi", "CHEBI:27902").name == "CHEBI_27902.json"
    assert cache.get("chebi", "CHEBI:27902") == {"name": "tetracycline"}
    assert cache.get("chebi", "CHEBI:1") is None


def test_get_or_fetch_only_fetches_once(tmp_path):
    cache = JsonFileCache(tmp_path)
    calls = []

    def fetch():
        calls.append(1)
        return {"v": 1}

    assert cache.get_or_fetch("ns", "k", fetch) == {"v": 1}
    assert cache.get_or_fetch("ns", "k", fetch) == {"v": 1}
    assert len(calls) == 1
    assert json.loads(cache.path("ns", "k").read_text())["v"] == 1


def test_offline_mode_raises_on_cache_miss(tmp_path, monkeypatch):
    cache = JsonFileCache(tmp_path)
    monkeypatch.setenv(OFFLINE_ENV, "1")
    with pytest.raises(NetworkUnavailable):
        cache.get_or_fetch("ns", "missing", lambda: {"never": True})
    cache.put("ns", "present", {"ok": True})
    assert cache.get_or_fetch("ns", "present", lambda: 1/0) == {"ok": True}
