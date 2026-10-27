"""Load ECOSTRESS spectral library files (header block + two-column data)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from leafrtm.db.adapters.core import DatasetTables, to_tables
from leafrtm.db.sources import ecostress
from leafrtm.db.sources.common import raw_dir


def parse(path):
    header, rows = {}, []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) == 2 and not rows and ":" not in line:
            rows.append(parts)
        elif rows:
            rows.append(parts)
        elif ":" in line:
            k, v = line.split(":", 1)
            header[k.strip()] = v.strip()
    data = np.array(rows, dtype=float)
    return header, data[:, 0], data[:, 1]


def load(cfg: dict) -> DatasetTables:
    files = sorted(p for c in ecostress.CLASSES for p in (raw_dir("ecostress") / c).glob("*.txt"))
    if not files:
        ecostress.fetch_all()
        files = sorted(p for c in ecostress.CLASSES for p in (raw_dir("ecostress") / c).glob("*.txt"))
    metas, spectra = [], []
    for f in files:
        h, x, y = parse(f)
        unit = h.get("X Units", "").lower()
        wl = x * 1000.0 if "micrometer" in unit else x
        h["file"] = f.name
        metas.append(h)
        spectra.append(pd.Series(y, index=wl).groupby(level=0).mean())
    wl = np.unique(np.concatenate([s.index.to_numpy() for s in spectra]))
    vals = np.full((len(spectra), len(wl)), np.nan)
    for i, s in enumerate(spectra):
        vals[i, np.searchsorted(wl, s.index.to_numpy())] = s.to_numpy()
    meta = pd.DataFrame(metas)
    origin = meta["Origin"].fillna("").str.split(";", expand=True)
    meta["lat"] = pd.to_numeric(origin[0], errors="coerce")
    meta["lon"] = pd.to_numeric(origin[1], errors="coerce") if origin.shape[1] > 1 else np.nan
    info = {"title": "ECOSTRESS spectral library v1.0", "url": ecostress.BASE,
            "doi": "10.1016/j.rse.2019.05.015",
            "citation": "Meerdink, S. K., Hook, S. J., Roberts, D. A., Abbott, E. A. (2019). The ECOSTRESS spectral library version 1.0. Remote Sens. Environ. 230, 111196.",
            "license": None}
    return to_tables(cfg, meta, wl, vals, info)
