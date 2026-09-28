"""Assembly-scoped entity extraction.

Which chemical entities are present in one biological assembly of a PDB entry:

* PDBe ``/pdb/entry/assembly/{pdb_id}`` lists, per assembly, the entities with their chains and
  copy numbers; PDBe ``/pdb/entry/molecules/{pdb_id}`` gives the CCD code of each bound entity.
* The entry's entities are queried in the PDBe search API (``/pdbe/search/pdb/select``), whose
  per-entity documents carry ``prd_id`` / ``prd_class`` / ``prd_name`` / ``prd_type`` for BIRD
  reference molecules.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import datetime as _dt

from .cache import HttpError, JsonFileCache

PDBE_ASSEMBLY_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/assembly/{pdb_id}"
PDBE_MOLECULES_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id}"
PDBE_SEARCH_URL = "https://www.ebi.ac.uk/pdbe/search/pdb/select"
SEARCH_ROWS = 500

KIND_CCD = "CCD"
KIND_PRD = "PRD"



def _now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat()


class EntryNotFound(LookupError):
    pass


class AssemblyNotFound(LookupError):
    pass


@dataclass
class AssemblyEntity:
    entity_id: str                      # CCD code (e.g. TAC) or PRD id (e.g. PRD_000226)
    entity_kind: str                    # "CCD" | "PRD"
    name: str | None = None             # deposited name
    pdb_entity_ids: list[int] = field(default_factory=list)
    chains: list[str] = field(default_factory=list)
    copies: int = 0

    def to_dict(self) -> dict[str, Any]:
        return {"entity_id": self.entity_id, "entity_kind": self.entity_kind, "name": self.name, "pdb_entity_ids": self.pdb_entity_ids, "chains": self.chains, "copies": self.copies}


class AssemblyClient:
    def __init__(self, transport, cache: JsonFileCache):
        self.transport = transport
        self.cache = cache

    # -- raw (cached) --------------------------------------------------------
    def _pdbe(self, url_tpl: str, namespace: str, pdb_id: str) -> list[dict[str, Any]]:
        pid = pdb_id.lower()

        def fetch():
            try:
                data = self.transport.get_json(url_tpl.format(pdb_id=pid))
            except HttpError as exc:
                if exc.status == 404:
                    return {}
                raise
            return data or {}

        data = self.cache.get_or_fetch(namespace, pid, fetch)
        if not data or pid not in data:
            raise EntryNotFound(f"PDB entry {pdb_id!r} not found at PDBe")
        return data[pid]

    def assemblies_raw(self, pdb_id: str) -> list[dict[str, Any]]:
        return self._pdbe(PDBE_ASSEMBLY_URL, "pdbe/assembly", pdb_id)

    def molecules_raw(self, pdb_id: str) -> list[dict[str, Any]]:
        return self._pdbe(PDBE_MOLECULES_URL, "pdbe/molecules", pdb_id)

    def prd_map(self, pdb_id: str) -> dict[int, dict[str, Any]]:
        """{entity_id: {"prd_id", "class", "name", "type"}} for every BIRD entity of the entry."""
        pid = pdb_id.lower()

        def fetch():
            docs: list[dict[str, Any]] = []
            start = 0
            while True:
                params = {"q": f"pdb_id:{pid}", "fl": "entity_id,prd_id,prd_class,prd_name,prd_type,molecule_type", "rows": SEARCH_ROWS, "start": start, "wt": "json"}
                resp = self.transport.get_json(PDBE_SEARCH_URL, params=params)
                page = (resp.get("response") or {}).get("docs") or []
                docs.extend(page)
                num_found = int((resp.get("response") or {}).get("numFound") or 0)
                start += SEARCH_ROWS
                if start >= num_found or not page:
                    break
            return {"docs": [d for d in docs if d.get("prd_id")], "retrieved_at": _now()}

        data = self.cache.get_or_fetch("pdbe/prd_entities", pid, fetch)
        out: dict[int, dict[str, Any]] = {}
        for d in data.get("docs") or []:
            prd = d.get("prd_id")
            prd = prd[0] if isinstance(prd, list) else prd
            if not prd or d.get("entity_id") is None:
                continue
            first = lambda v: v[0] if isinstance(v, list) and v else (v if not isinstance(v, list) else None)  # noqa: E731
            out[int(d["entity_id"])] = {"prd_id": prd, "class": first(d.get("prd_class")), "name": first(d.get("prd_name")), "type": first(d.get("prd_type"))}
        return out

    # -- public --------------------------------------------------------------
    def assembly_ids(self, pdb_id: str) -> list[str]:
        return [str(a["assembly_id"]) for a in self.assemblies_raw(pdb_id)]

    def entities(self, pdb_id: str, assembly_id: str | int) -> list[AssemblyEntity]:
        """CCD and PRD entities present in the given assembly (water excluded)."""
        aid = str(assembly_id).strip()
        assemblies = self.assemblies_raw(pdb_id)
        match = [a for a in assemblies if str(a.get("assembly_id")) == aid]
        if not match:
            raise AssemblyNotFound(f"assembly {aid!r} not found for {pdb_id.upper()}; available: {', '.join(str(a['assembly_id']) for a in assemblies)}")
        asm = match[0]
        molecules = {int(m["entity_id"]): m for m in self.molecules_raw(pdb_id)}

        ccd: dict[str, AssemblyEntity] = {}
        polymer_rows: dict[int, dict[str, Any]] = {}
        for e in asm.get("entities") or []:
            eid = int(e["entity_id"])
            mtype = (e.get("molecule_type") or "").lower()
            if mtype == "water":
                continue
            if mtype == "bound":
                mol = molecules.get(eid, {})
                for comp in mol.get("chem_comp_ids") or []:
                    ent = ccd.setdefault(comp.upper(), AssemblyEntity(entity_id=comp.upper(), entity_kind=KIND_CCD, name=(mol.get("molecule_name") or [None])[0]))
                    ent.pdb_entity_ids.append(eid)
                    ent.chains.extend(e.get("in_chains") or [])
                    ent.copies += int(e.get("number_of_copies") or len(e.get("in_chains") or []))
                continue
            polymer_rows[eid] = e

        prd: dict[str, AssemblyEntity] = {}
        for eid, info in self.prd_map(pdb_id).items():
            e = polymer_rows.get(eid)
            if e is None:
                continue  # BIRD entity exists in the entry but not in this assembly
            ent = prd.setdefault(info["prd_id"], AssemblyEntity(entity_id=info["prd_id"], entity_kind=KIND_PRD, name=info.get("name") or (e.get("molecule_name") or [None])[0]))
            ent.pdb_entity_ids.append(eid)
            ent.chains.extend(e.get("in_chains") or [])
            ent.copies += int(e.get("number_of_copies") or len(e.get("in_chains") or []))
        return sorted(ccd.values(), key=lambda x: x.entity_id) + sorted(prd.values(), key=lambda x: x.entity_id)
