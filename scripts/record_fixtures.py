"""Record real API responses as replayable test fixtures (VCR style).

Run once (network required):  python scripts/record_fixtures.py
Writes tests/fixtures/http/index.json + one JSON file per request. Tests replay these through
``FakeTransport`` and never touch the network.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibiotic_annotation.cache import HttpTransport, JsonFileCache  # noqa: E402
from antibiotic_annotation.chebi import ChebiClient, ols4_term_url, trim_backend, trim_ols4  # noqa: E402
from antibiotic_annotation.identity import RcsbClient  # noqa: E402
from antibiotic_annotation.mapping import Mapper, UniChemClient  # noqa: E402

FIX = ROOT / "tests" / "fixtures" / "http"

IDS = ["PRD_000505", "PRD_000193", "5I0", "GLC", "TAC", "KSG", "SCM", "GTP", "ORN", "TPF", "X8Q", "FYG", "T1C", "HY0", "6UQ", "NMY", "LLL", "VIR", "KIR", "PRD_000226", "EM1", "HOH", "CHEBI:26710", "ZZZZ9"]
ASSEMBLIES = [("5J7L", 1), ("5J7L", 2), ("4V7T", 1), ("4V7T", 2), ("4V7U", 1), ("4V85", 1), ("4U1U", 1), ("4V64", 1), ("5AFI", 1), ("1ATP", 1)]
EXTRA_CHEBI = ["CHEBI:36047", "CHEBI:33282", "CHEBI:33281", "CHEBI:22507", "CHEBI:25105", "CHEBI:27933", "CHEBI:26895", "CHEBI:7507", "CHEBI:17833", "CHEBI:87209", "CHEBI:48947", "CHEBI:48001", "CHEBI:24835", "CHEBI:9999999"]


def request_key(method: str, url: str, payload) -> str:
    return hashlib.sha1(json.dumps([method, url, payload], sort_keys=True).encode()).hexdigest()[:16]


def bundle_name(url: str) -> str:
    if "pdbe" in url:
        return "pdbe"
    if "ols4" in url:
        return "ols4"
    if "chebi/backend" in url:
        return "chebi_backend"
    if "rcsb" in url:
        return "rcsb"
    return "unichem"


class RecordingTransport(HttpTransport):
    def __init__(self):
        super().__init__()
        FIX.mkdir(parents=True, exist_ok=True)
        self.bundles: dict[str, dict] = {}

    def _bundle(self, name):
        if name not in self.bundles:
            path = FIX / f"{name}.json"
            self.bundles[name] = json.loads(path.read_text()) if path.exists() else {}
        return self.bundles[name]

    def flush(self):
        for name, b in self.bundles.items():
            (FIX / f"{name}.json").write_text(json.dumps(b, indent=0, sort_keys=True, ensure_ascii=False))

    def _record(self, method, url, payload, fn):
        key = request_key(method, url, payload)
        from antibiotic_annotation.cache import HttpError

        try:
            data = fn()
            status = 200
        except HttpError as exc:
            data, status = {"error": str(exc)}, exc.status
        # trim the big ontology payloads
        if status == 200 and "ols4/api/v2/ontologies/chebi/classes" in url:
            data = trim_ols4(data)
        elif status == 200 and "chebi/backend/api/public/compound" in url:
            data = trim_backend(data)
        self._bundle(bundle_name(url))[key] = {"method": method, "url": url, "payload": payload, "status": status, "data": data}
        if status != 200:
            raise HttpError(url, status)
        return data

    def get_json(self, url, params=None):
        return self._record("GET", url, params, lambda: super(RecordingTransport, self).get_json(url, params))

    def post_json(self, url, payload):
        return self._record("POST", url, payload, lambda: super(RecordingTransport, self).post_json(url, payload))


def main():
    import tempfile

    from antibiotic_annotation.identity import IdentityNotFound

    tr = RecordingTransport()
    cache = JsonFileCache(tempfile.mkdtemp())
    rcsb, uni, chebi = RcsbClient(tr, cache), UniChemClient(tr, cache), ChebiClient(tr, cache)
    from antibiotic_annotation.pipeline import Pipeline
    pipe = Pipeline(cache_dir=cache.root, transport=tr, synonyms=None, related_parents=None)
    mapper = pipe.mapper
    chebi.release()
    to_expand: set[str] = set(EXTRA_CHEBI)
    for i in IDS:
        try:
            ident = rcsb.identity(i)
        except IdentityNotFound:
            print(i, "-> not found (recorded)")
            continue
        m = mapper.map(ident)
        print(i, ident.entity_kind, ident.inchikey, m.status, m.method, m.primary_chebi_id, m.equivalent_chebi_ids)
        to_expand.update(m.chebi_ids_for_evidence)
        to_expand.update(c.chebi_id for c in m.candidates)
    for cid in sorted(to_expand):
        if not chebi.has_term(cid):
            print(cid, "missing")
            continue
        roles = chebi.inherited_roles(cid)
        chebi.role_closure(roles.keys())
        for r in chebi.term(cid).relations + chebi.term(cid).incoming:
            chebi.term(r.target_id)
        print(cid, chebi.name(cid), "ancestors", len(chebi.ancestors(cid)), "roles", len(roles))
    tr.flush()
    from antibiotic_annotation.api import inspect_assembly
    from antibiotic_annotation.assembly import AssemblyClient, AssemblyNotFound, EntryNotFound

    for pdb, asm in ASSEMBLIES:
        r = inspect_assembly(pdb, asm, pipeline=pipe)
        print(pdb, asm, "inspected", [e.entity_id for e in r.entities_inspected], "hits", [h.entity_id for h in r.hits])
    for pdb, asm in [("5J7L", 9), ("XXXX", 1)]:
        try:
            AssemblyClient(tr, cache).entities(pdb, asm)
        except (AssemblyNotFound, EntryNotFound) as exc:
            print("expected error recorded:", exc)
    tr.flush()
    print("fixtures:", {k: len(v) for k, v in tr.bundles.items()})


if __name__ == "__main__":
    main()
