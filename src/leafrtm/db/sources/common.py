"""Shared download helpers."""

from __future__ import annotations

import gzip
import os
import shutil
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[4]


def data_dir() -> Path:
    return Path(os.environ.get("LEAFRTM_DATA", ROOT / "data"))


def raw_dir(source: str) -> Path:
    d = data_dir() / "raw" / source
    d.mkdir(parents=True, exist_ok=True)
    return d


def get(url, *, session=None, retries=6, timeout=120, **kwargs):
    """GET with exponential backoff on network errors and 5xx."""
    s = session or requests
    for attempt in range(retries + 1):
        try:
            r = s.get(url, timeout=timeout, **kwargs)
            # 403/429 from EcoSIS are transient rate limits; retry those too
            if r.status_code < 500 and r.status_code not in (403, 429):
                r.raise_for_status()
                return r
        except (requests.ConnectionError, requests.Timeout):
            if attempt == retries:
                raise
        time.sleep(3 * 2**attempt)
    r.raise_for_status()
    return r


def download(url, dest: Path, *, session=None, gz=False, **kwargs) -> Path:
    """Stream ``url`` to ``dest`` (optionally gzip-compressed); atomic and cached."""
    dest = Path(dest)
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    s = session or requests
    for attempt in range(5):
        try:
            with s.get(url, stream=True, timeout=300, **kwargs) as r:
                r.raise_for_status()
                opener = gzip.open if gz else open
                with opener(tmp, "wb") as f:
                    shutil.copyfileobj(r.raw if not r.headers.get("content-encoding") else
                                       _iter_raw(r), f)
            break
        except (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError):
            if attempt == 4:
                raise
            time.sleep(2 ** (attempt + 1))
    tmp.rename(dest)
    return dest


class _iter_raw:
    """File-like wrapper decoding content-encoding transparently."""

    def __init__(self, r):
        self._it = r.iter_content(1 << 20)
        self._buf = b""

    def read(self, n=-1):
        while n < 0 or len(self._buf) < n:
            try:
                self._buf += next(self._it)
            except StopIteration:
                break
        if n < 0:
            out, self._buf = self._buf, b""
        else:
            out, self._buf = self._buf[:n], self._buf[n:]
        return out
