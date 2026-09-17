"""Phase 3, step 1: PDB chemical identity from the RCSB Data API.

Both Chemical Component Dictionary ids (``TAC``) and BIRD reference-molecule ids
(``PRD_000226``) are served by the same ``core/chemcomp`` endpoint. For BIRD entries the
``pdbx_reference_molecule`` block (class, type) is kept as an independent evidence channel.
"""
from __future__ import annotations

import datetime as _dt
import re
from typing import Any

from .cache import HttpError, JsonFileCache, NetworkUnavailable
from .models import ENTITY_BIRD, ENTITY_CCD, ENTITY_CHEBI_INPUT, ChemicalIdentity, normalise_chebi_id

RCSB_CHEMCOMP_URL = "https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}"
CACHE_NS = "rcsb/chemcomp"

_PRD_RE = re.compile(r"^PRD_\d{6}$")
_CCD_RE = re.compile(r"^[A-Za-z0-9]{1,5}$")


class IdentityNotFound(LookupError):
    """The id is not a known CCD / PRD entry."""


def classify_input_id(input_id: str) -> str:
    s = input_id.strip()
    if normalise_chebi_id(s) and s.upper().startswith("CHEBI"):
        return ENTITY_CHEBI_INPUT
    if _PRD_RE.match(s.upper()):
        return ENTITY_BIRD
    if _CCD_RE.match(s):
        return ENTITY_CCD
    raise ValueError(f"unrecognised input id: {input_id!r}")


def _now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat()


def parse_chemcomp(raw: dict[str, Any], input_id: str) -> ChemicalIdentity:
    cc = raw.get("chem_comp") or {}
    desc = raw.get("rcsb_chem_comp_descriptor") or {}
    ident = ChemicalIdentity(input_id=input_id)
    ident.name = cc.get("name")
    ident.ccd_type = cc.get("type")
    ident.formula = cc.get("formula")
    charge = cc.get("pdbx_formal_charge")
    ident.formal_charge = int(charge) if charge is not None else None
    ident.inchi = desc.get("InChI")
    ident.inchikey = desc.get("InChIKey")
    ident.smiles = desc.get("SMILES_stereo") or desc.get("SMILES")
    syn: list[str] = []
    for s in raw.get("rcsb_chem_comp_synonyms") or []:
        n = s.get("name")
        if n and n not in syn:
            syn.append(n)
    ident.synonyms = syn
    xrefs: dict[str, list[str]] = {}
    for r in raw.get("rcsb_chem_comp_related") or []:
        res, acc = r.get("resource_name"), r.get("resource_accession_code")
        if res and acc:
            xrefs.setdefault(res, [])
            if acc not in xrefs[res]:
                xrefs[res].append(acc)
    ident.xrefs = xrefs
    ident.atc_codes = sorted({a["annotation_id"] for a in raw.get("rcsb_chem_comp_annotation") or [] if a.get("type") == "ATC" and a.get("annotation_id")})
    prd = raw.get("pdbx_reference_molecule")
    if prd or (cc.get("id", "").upper().startswith("PRD_")):
        ident.entity_kind = ENTITY_BIRD
        if prd:
            ident.bird_class = prd.get("class")
            ident.bird_type = prd.get("type")
            ident.name = prd.get("name") or ident.name
    else:
        ident.entity_kind = ENTITY_CCD
    ident.retrieved_at = raw.get("_retrieved_at") or _now()
    return ident


class RcsbClient:
    def __init__(self, transport, cache: JsonFileCache):
        self.transport = transport
        self.cache = cache

    def raw_chemcomp(self, comp_id: str) -> dict[str, Any]:
        comp_id = comp_id.strip().upper()

        def fetch():
            try:
                data = self.transport.get_json(RCSB_CHEMCOMP_URL.format(comp_id=comp_id))
            except HttpError as exc:
                if exc.status == 404:
                    raise IdentityNotFound(comp_id) from exc
                raise
            data["_retrieved_at"] = _now()
            return data

        return self.cache.get_or_fetch(CACHE_NS, comp_id, fetch)

    def identity(self, input_id: str) -> ChemicalIdentity:
        kind = classify_input_id(input_id)
        if kind == ENTITY_CHEBI_INPUT:
            return ChemicalIdentity(input_id=normalise_chebi_id(input_id), entity_kind=ENTITY_CHEBI_INPUT, source="input", retrieved_at=_now())
        raw = self.raw_chemcomp(input_id)
        return parse_chemcomp(raw, input_id.strip().upper())
