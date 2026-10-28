"""EcoSIS (ecosis.org) leaf spectra packages.

Packages are discovered through the search API (Target Type = leaf) and
exported as wide CSV (metadata columns, then one column per wavelength).
"""

from __future__ import annotations

import json

import requests

from leafrtm.db.sources.common import download, get, raw_dir

API = "https://ecosis.org/api"


def search_leaf_packages(page=100) -> list[dict]:
    """All EcoSIS packages with Target Type = leaf (package JSON, abridged)."""
    out, start = [], 0
    while True:
        r = get(f"{API}/package/search", params={
            "text": "", "filters": json.dumps([{"Target Type": "leaf"}]),
            "start": start, "stop": start + page})
        d = r.json()
        out.extend(d["items"])
        start += page
        if start >= d["total"]:
            return out


def package_metadata(name: str, session=None) -> dict:
    path = raw_dir("ecosis") / f"{name}.json"
    if path.exists():
        return json.loads(path.read_text())
    d = get(f"{API}/package/{name}", session=session).json()
    path.write_text(json.dumps(d, indent=1))
    return d


def export(name: str, session=None):
    """Download a package's spectra + metadata export (gzipped CSV)."""
    return download(f"{API}/package/{name}/export?metadata=true",
                    raw_dir("ecosis") / f"{name}.csv.gz", session=session, gz=True)


def fetch_all(names) -> None:
    s = requests.Session()
    for n in names:
        package_metadata(n, s)
        export(n, s)
