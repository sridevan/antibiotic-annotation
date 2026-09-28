"""Shared fixtures: a replaying HTTP transport fed from recorded real responses.

No test touches the network. ``FakeTransport`` raises ``NetworkUnavailable`` for any request
that was not recorded, which doubles as a network-failure simulation.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from antibiotic_annotation.cache import HttpError, JsonFileCache, NetworkUnavailable
from antibiotic_annotation.chebi import ChebiClient
from antibiotic_annotation.identity import CompoundClient
from antibiotic_annotation.mapping import ChebiLookup, Mapper, UniChemClient

FIXTURES = Path(__file__).parent / "fixtures" / "http"


def _load_bundles() -> dict:
    out: dict = {}
    for f in sorted(FIXTURES.glob("*.json")):
        out.update(json.loads(f.read_text(encoding="utf-8")))
    return out


_BUNDLES = _load_bundles()


def request_key(method: str, url: str, payload) -> str:
    return hashlib.sha1(json.dumps([method, url, payload], sort_keys=True).encode()).hexdigest()[:16]


class FakeTransport:
    def __init__(self, fail_all: bool = False):
        self.calls: list[tuple[str, str]] = []
        self.fail_all = fail_all

    def _do(self, method, url, payload):
        self.calls.append((method, url))
        if self.fail_all:
            raise NetworkUnavailable(f"simulated network failure for {url}")
        rec = _BUNDLES.get(request_key(method, url, payload))
        if rec is None:
            raise NetworkUnavailable(f"no recorded fixture for {method} {url} {payload}")
        if rec["status"] != 200:
            raise HttpError(url, rec["status"])
        return json.loads(json.dumps(rec["data"]))  # defensive copy

    def get_json(self, url, params=None):
        return self._do("GET", url, params)

    def post_json(self, url, payload):
        return self._do("POST", url, payload)


class FailingTransport(FakeTransport):
    def __init__(self):
        super().__init__(fail_all=True)


@pytest.fixture
def transport():
    return FakeTransport()


@pytest.fixture
def cache(tmp_path):
    return JsonFileCache(tmp_path / "cache")


@pytest.fixture
def compounds(transport, cache):
    return CompoundClient(transport, cache)


@pytest.fixture
def chebi(transport, cache):
    return ChebiClient(transport, cache)


@pytest.fixture
def mapper(transport, cache, chebi):
    def lookup(c):
        t = chebi.term(c)
        forms = {r.target_id for r in t.relations + t.incoming if r.relation in ("is_conjugate_acid_of", "is_conjugate_base_of", "is_tautomer_of")}
        return ChebiLookup(stars=t.stars, inchikey=t.inchikey, is_class_like=False, num_descendants=t.num_descendants, forms=frozenset(forms))

    return Mapper(UniChemClient(transport, cache), lookup=lookup)
