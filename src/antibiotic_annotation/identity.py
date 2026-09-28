"""Deposited chemical identity from the PDBe compound API.

``GET https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{id}`` serves both Chemical Component
Dictionary codes (``TAC``) and BIRD reference-molecule ids (``PRD_000226``). For BIRD entries the
``compound_type`` ("Oligopeptide", ...) and ``compound_classes`` (e.g. ["antibiotic"]) are kept
as an independent evidence channel.
"""
from __future__ import annotations

import datetime as _dt
import re
from typing import Any

from .cache import HttpError, JsonFileCache
from .models import ENTITY_BIRD, ENTITY_CCD, ENTITY_CHEBI_INPUT, ChemicalIdentity, normalise_chebi_id

PDBE_COMPOUND_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{comp_id}"
CACHE_NS = "pdbe/compound"

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


def _pick_smiles(entries: list[dict[str, Any]] | None) -> str | None:
    """Prefer a stereo-aware SMILES (contains @) from any program, else the first one."""
    entries = entries or []
    for e in entries:
        if "@" in (e.get("name") or ""):
            return e["name"]
    return entries[0].get("name") if entries else None


def parse_compound(raw: dict[str, Any], input_id: str) -> ChemicalIdentity:
    ident = ChemicalIdentity(input_id=input_id, source="PDBe")
    ident.name = raw.get("name")
    ident.ccd_type = raw.get("compound_type")
    ident.formula = raw.get("formula")
    charge = raw.get("formal_charge")
    try:
        ident.formal_charge = int(charge) if charge not in (None, "") else None
    except (TypeError, ValueError):
        ident.formal_charge = None
    ident.inchi = raw.get("inchi")
    ident.inchikey = raw.get("inchi_key")
    ident.smiles = _pick_smiles(raw.get("smiles"))
    syn: list[str] = []
    for s in raw.get("synonyms") or []:
        v = s.get("value") if isinstance(s, dict) else s
        if v and v not in syn:
            syn.append(v)
    ident.synonyms = syn
    xrefs: dict[str, list[str]] = {}
    for x in raw.get("cross_links") or []:
        res, acc = x.get("resource"), x.get("resource_id")
        if res and acc is not None:
            xrefs.setdefault(res, [])
            if str(acc) not in xrefs[res]:
                xrefs[res].append(str(acc))
    ident.xrefs = xrefs
    if input_id.upper().startswith("PRD_"):
        ident.entity_kind = ENTITY_BIRD
        classes = raw.get("compound_classes") or []
        ident.bird_class = ", ".join(str(c) for c in classes) if classes else None
        ident.bird_type = raw.get("compound_type")
    else:
        ident.entity_kind = ENTITY_CCD
    ident.retrieved_at = raw.get("_retrieved_at") or _now()
    return ident


class CompoundClient:
    def __init__(self, transport, cache: JsonFileCache):
        self.transport = transport
        self.cache = cache

    def raw_compound(self, comp_id: str) -> dict[str, Any]:
        comp_id = comp_id.strip().upper()

        def fetch():
            try:
                data = self.transport.get_json(PDBE_COMPOUND_URL.format(comp_id=comp_id))
            except HttpError as exc:
                if exc.status == 404:
                    data = {}
                else:
                    raise
            entries = (data or {}).get(comp_id) or []
            if not entries:
                raise IdentityNotFound(comp_id)
            rec = dict(entries[0])
            rec["_retrieved_at"] = _now()
            return rec

        return self.cache.get_or_fetch(CACHE_NS, comp_id, fetch)

    def identity(self, input_id: str) -> ChemicalIdentity:
        kind = classify_input_id(input_id)
        if kind == ENTITY_CHEBI_INPUT:
            return ChemicalIdentity(input_id=normalise_chebi_id(input_id), entity_kind=ENTITY_CHEBI_INPUT, source="input", retrieved_at=_now())
        raw = self.raw_compound(input_id)
        return parse_compound(raw, input_id.strip().upper())
