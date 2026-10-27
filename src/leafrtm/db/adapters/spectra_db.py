"""Load a spectra_db dataset directory (``metadata.csvy`` + ``spectra/*.csvy``).

Each spectra file holds one observation: a ``wavelengths`` column and one
column per spectrum; its YAML header lists ``spectra_types`` (one value, or
one per column). The loader returns one metadata row per spectrum, with the
observation's metadata plus ``spectra_type`` and ``spectrum`` (column name).
"""

from __future__ import annotations

import io

import numpy as np
import pandas as pd
import yaml

from leafrtm.db.adapters.core import DatasetTables, to_tables
from leafrtm.db.sources import spectra_db


def read_csvy(path):
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    header = [ln[1:] for ln in lines if ln.startswith("#")]
    body = "\n".join(ln for ln in lines if not ln.startswith("#"))
    meta = yaml.safe_load("\n".join(h for h in header if h.strip() != "---")) or {}
    return meta, pd.read_csv(io.StringIO(body), low_memory=False)


def load_dir(name):
    root = spectra_db.checkout() / name
    header, meta = read_csvy(root / "metadata.csvy")
    meta = meta.set_index("observation_id", drop=False)
    rows, specs = [], []
    for f in sorted((root / "spectra").glob("*.csvy")):
        h, df = read_csvy(f)
        df = df.dropna(subset=["wavelengths"]).sort_values("wavelengths")
        cols = [c for c in df.columns if c != "wavelengths"]
        types = h.get("spectra_types", "R")
        types = types if isinstance(types, list) else [types] * len(cols)
        obs = f.stem
        base = meta.loc[obs].to_dict() if obs in meta.index else {"observation_id": obs}
        for c, t in zip(cols, types):
            rows.append({**base, "spectra_type": t, "spectrum": c})
            specs.append(pd.Series(pd.to_numeric(df[c], errors="coerce").to_numpy(),
                                   index=df["wavelengths"].to_numpy(float)))
    wl = np.unique(np.concatenate([s.index.to_numpy() for s in specs]))
    vals = np.full((len(specs), len(wl)), np.nan)
    for i, s in enumerate(specs):
        s = s.groupby(level=0).mean()
        vals[i, np.searchsorted(wl, s.index.to_numpy())] = s.to_numpy()
    return header, pd.DataFrame(rows), wl, vals


def load(cfg: dict) -> DatasetTables:
    header, meta, wl, vals = load_dir(cfg["source_id"])
    info = {"title": header.get("long_name"), "url": header.get("URL"),
            "doi": header.get("doi"), "citation": None, "license": None}
    info.update({k: cfg[k] for k in ("title", "url", "doi", "citation") if k in cfg})
    return to_tables(cfg, meta, wl, vals, info)
