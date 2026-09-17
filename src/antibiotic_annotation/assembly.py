"""Assembly-scoped entity extraction.

Which chemical entities are present in one biological assembly of a PDB entry:

* PDBe ``/pdb/entry/assembly/{pdb_id}`` lists, per assembly, the entities with their chains and
  copy numbers; PDBe ``/pdb/entry/molecules/{pdb_id}`` gives the CCD code of each bound entity.
* PDBe does not expose BIRD/PRD ids, so the (short) polymer entities of the assembly are looked
  up in RCSB GraphQL ``polymer_entities`` for ``prd_id``. PDBe and RCSB share mmCIF entity ids.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .cache import HttpError, JsonFileCache

PDBE_ASSEMBLY_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/assembly/{pdb_id}"
PDBE_MOLECULES_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id}"
RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"

KIND_CCD = "CCD"
KIND_PRD = "PRD"

# Polymer entities longer than this are never BIRD peptide-like molecules; skipping them keeps
# the RCSB lookup small for ribosomes (which have ~50 protein/RNA entities).
MAX_BIRD_POLYMER_LENGTH = 100


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

    def prd_map(self, pdb_id: str, polymer_entity_ids: list[int]) -> dict[int, dict[str, Any]]:
        """{entity_id: {"prd_id": ..., "description": ...}} for the polymer entities that carry a PRD id."""
        pid = pdb_id.upper()
        ids = sorted(set(polymer_entity_ids))
        if not ids:
            return {}
        key = f"{pid}_{'_'.join(str(i) for i in ids)}"
        query = ("{ polymer_entities(entity_ids:[%s]) { rcsb_id rcsb_polymer_entity { pdbx_description } "
                 "rcsb_polymer_entity_container_identifiers { entity_id prd_id } } }") % ",".join(f'"{pid}_{i}"' for i in ids)
        data = self.cache.get_or_fetch("rcsb/prd_map", key, lambda: self.transport.post_json(RCSB_GRAPHQL_URL, {"query": query}))
        out: dict[int, dict[str, Any]] = {}
        for pe in (data.get("data") or {}).get("polymer_entities") or []:
            if not pe:
                continue
            ci = pe.get("rcsb_polymer_entity_container_identifiers") or {}
            if ci.get("prd_id"):
                out[int(ci["entity_id"])] = {"prd_id": ci["prd_id"], "description": (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description")}
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
        polymer_ids: list[int] = []
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
            length = molecules.get(eid, {}).get("length")
            if length is not None and length > MAX_BIRD_POLYMER_LENGTH:
                continue
            polymer_ids.append(eid)
            polymer_rows[eid] = e

        prd: dict[str, AssemblyEntity] = {}
        for eid, info in self.prd_map(pdb_id, polymer_ids).items():
            e = polymer_rows[eid]
            ent = prd.setdefault(info["prd_id"], AssemblyEntity(entity_id=info["prd_id"], entity_kind=KIND_PRD, name=info.get("description")))
            ent.pdb_entity_ids.append(eid)
            ent.chains.extend(e.get("in_chains") or [])
            ent.copies += int(e.get("number_of_copies") or len(e.get("in_chains") or []))
        return sorted(ccd.values(), key=lambda x: x.entity_id) + sorted(prd.values(), key=lambda x: x.entity_id)
