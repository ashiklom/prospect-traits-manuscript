"""Load EcoSIS package exports for a dataset config."""

from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd

from leafrtm.db.adapters.core import DatasetTables, to_tables
from leafrtm.db.sources import ecosis
from leafrtm.db.sources.common import raw_dir


def _wl(c) -> float | None:
    """Wavelength from a column name ("512", "512.5", "512 nm"), else None."""
    try:
        return float(str(c).strip().removesuffix("nm").strip())
    except ValueError:
        return None


def _is_wl(c) -> bool:
    return _wl(c) is not None


_BAND = re.compile(r"^B(\d{3,4}(?:\.\d+)?)$")


def read_export(package: str) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Metadata, wavelengths and values of an export.

    Some exports carry part of the spectrum as metadata-like columns named
    ``B<nm>``; those are read as bands too (and also kept in the metadata,
    where configs may use them, e.g. as spectral fingerprints).
    """
    path = raw_dir("ecosis") / f"{package}.csv.gz"
    if not path.exists():
        ecosis.export(package)
    df = pd.read_csv(path, low_memory=False, encoding="utf-8-sig")
    wl_cols = [c for c in df.columns if _is_wl(c)]
    b_cols = {c: float(_BAND.match(c.strip()).group(1)) for c in df.columns
              if isinstance(c, str) and _BAND.match(c.strip())}
    b_cols = {c: w for c, w in b_cols.items() if 300 <= w <= 2600
              and w not in {_wl(x) for x in wl_cols}}
    meta = df.drop(columns=wl_cols)
    meta.columns = [c.strip() for c in meta.columns]
    cols = wl_cols + list(b_cols)
    vals = df[cols].apply(pd.to_numeric, errors="coerce").to_numpy(np.float64)
    wl = np.array([_wl(c) for c in wl_cols] + list(b_cols.values()))
    return meta, wl, vals


def package_info(package: str) -> dict:
    m = ecosis.package_metadata(package)
    e = m.get("ecosis", {})
    cit = m.get("Citation")
    doi = m.get("Citation DOI") or e.get("doi")
    return {
        "title": e.get("package_title"),
        "url": f"https://ecosis.org/package/{package}",
        "citation": "; ".join(cit) if isinstance(cit, list) else cit,
        "doi": "; ".join(doi) if isinstance(doi, list) else doi,
        "license": e.get("license"),
    }


def load(cfg: dict) -> DatasetTables:
    packages = cfg["source_id"] if isinstance(cfg["source_id"], list) else [cfg["source_id"]]
    metas, wls, vals = [], [], []
    drop_int = cfg.get("spectra", {}).get("drop_integer_bands", False)
    for i, p in enumerate(packages):
        m, w, v = read_export(p)
        if drop_int is True or (isinstance(drop_int, list) and i in drop_int):
            keep = w != np.round(w)
            w, v = w[keep], v[:, keep]
        m["_package"] = p
        m["_package_index"] = i
        metas.append(m)
        wls.append(w)
        vals.append(v)
    if len(packages) == 1:
        meta, wl, val = metas[0], wls[0], vals[0]
    else:  # union of wavelength grids
        wl = np.unique(np.concatenate(wls))
        val = np.full((sum(len(v) for v in vals), len(wl)), np.nan)
        r = 0
        for w, v in zip(wls, vals):
            val[r:r + len(v), np.searchsorted(wl, w)] = v
            r += len(v)
        meta = pd.concat(metas, ignore_index=True)
    info = package_info(packages[0])
    info["citation"] = cfg.get("citation") or info["citation"]
    if len(packages) > 1:
        info["url"] = "; ".join(f"https://ecosis.org/package/{p}" for p in packages)
    return to_tables(cfg, meta, wl, val, info)


def raw_metadata_json(package: str) -> str:
    return json.dumps(ecosis.package_metadata(package), indent=1)
