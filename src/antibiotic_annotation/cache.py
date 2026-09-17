"""JSON file cache and a thin HTTP transport with retries and an offline switch.

All network access in the package goes through :class:`HttpTransport`. Tests inject a
fake transport; the offline switch (env ``ANTIBIOTIC_ANNOTATION_OFFLINE=1``) makes any
cache miss raise :class:`NetworkUnavailable` instead of touching the network.
"""
from __future__ import annotations

import json
import os
import re
import time
from pathlib import Path
from typing import Any, Callable

OFFLINE_ENV = "ANTIBIOTIC_ANNOTATION_OFFLINE"


class NetworkUnavailable(RuntimeError):
    """Raised when a resource is not cached and the network cannot (or must not) be used."""


class HttpError(RuntimeError):
    def __init__(self, url: str, status: int | None, message: str = ""):
        super().__init__(f"HTTP {status} for {url} {message}".strip())
        self.url = url
        self.status = status


def is_offline() -> bool:
    return os.environ.get(OFFLINE_ENV, "").strip() not in ("", "0", "false", "False")


def _safe_key(key: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.:-]", "_", key).replace(":", "_")


class JsonFileCache:
    """One JSON file per (namespace, key) under ``root``."""

    def __init__(self, root: str | os.PathLike):
        self.root = Path(root)

    def path(self, namespace: str, key: str) -> Path:
        return self.root / namespace / f"{_safe_key(key)}.json"

    def get(self, namespace: str, key: str) -> Any | None:
        p = self.path(namespace, key)
        if not p.exists():
            return None
        with p.open("r", encoding="utf-8") as fh:
            return json.load(fh)

    def put(self, namespace: str, key: str, value: Any) -> None:
        p = self.path(namespace, key)
        p.parent.mkdir(parents=True, exist_ok=True)
        tmp = p.with_suffix(".json.tmp")
        with tmp.open("w", encoding="utf-8") as fh:
            json.dump(value, fh, indent=1, sort_keys=True, ensure_ascii=False)
        tmp.replace(p)

    def has(self, namespace: str, key: str) -> bool:
        return self.path(namespace, key).exists()

    def get_or_fetch(self, namespace: str, key: str, fetch: Callable[[], Any]) -> Any:
        """Return the cached value; otherwise call ``fetch`` (unless offline) and cache it."""
        cached = self.get(namespace, key)
        if cached is not None:
            return cached
        if is_offline():
            raise NetworkUnavailable(f"{namespace}/{key} not cached and offline mode is on")
        value = fetch()
        self.put(namespace, key, value)
        return value


class HttpTransport:
    """Small wrapper around ``requests`` with retries. Tests replace it with a fake."""

    def __init__(self, timeout: float = 30.0, retries: int = 3, backoff: float = 1.0, user_agent: str = "antibiotic-annotation/0.1"):
        self.timeout = timeout
        self.retries = retries
        self.backoff = backoff
        self.user_agent = user_agent
        self._session = None

    def _sess(self):
        if self._session is None:
            import requests  # imported lazily so tests never need the network stack

            self._session = requests.Session()
            self._session.headers["User-Agent"] = self.user_agent
            self._session.headers["Accept"] = "application/json"
        return self._session

    def _request(self, method: str, url: str, **kw) -> Any:
        if is_offline():
            raise NetworkUnavailable(f"offline mode: refusing {method} {url}")
        import requests

        last_exc: Exception | None = None
        for attempt in range(self.retries):
            try:
                resp = self._sess().request(method, url, timeout=self.timeout, **kw)
            except requests.RequestException as exc:  # connection errors, timeouts
                last_exc = exc
                time.sleep(self.backoff * (attempt + 1))
                continue
            if resp.status_code == 404:
                raise HttpError(url, 404, "not found")
            if resp.status_code >= 500 or resp.status_code == 429:
                last_exc = HttpError(url, resp.status_code)
                time.sleep(self.backoff * (attempt + 1))
                continue
            if resp.status_code >= 400:
                raise HttpError(url, resp.status_code, resp.text[:200])
            try:
                return resp.json()
            except ValueError as exc:
                raise HttpError(url, resp.status_code, "non-JSON response") from exc
        raise NetworkUnavailable(f"{method} {url} failed after {self.retries} attempts: {last_exc}")

    def get_json(self, url: str, params: dict[str, Any] | None = None) -> Any:
        return self._request("GET", url, params=params)

    def post_json(self, url: str, payload: Any) -> Any:
        return self._request("POST", url, json=payload)
