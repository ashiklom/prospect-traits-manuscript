"""Load the Yang et al. (2016) leaf spectra and traits from the authors' files.

Layout of the raw files (``sources/yang2016.py``):

* ``<year>_<site>_leaftraits.xlsx``: one sheet per trait (Chla, Chlb, TotChl, Car,
  LMA, TotC, TotN; capitalization differs between files). Rows are sampling days
  (first column = day of year), columns are leaf samples. Most sheets have a
  header row of sample names; the 2011 TotC sheet has none and the 2011 TotN
  header gives only the sample group (``AUN`` for AUN1-3), so sample names come
  from the Chla sheet and the other sheets are matched by column position
  (headers, where present, are checked against them).
* Spectra: one CSV per sampling day, no header; column 1 is wavelength (trailing
  rows with wavelength 0 are padding), every other column is one leaf sample.
  2011 Martha's Vineyard: ``2011-MV/NNsmoothed.csv`` (R only); file NN is the
  NN-th sampling day of the traits workbook. 2012 Harvard Forest:
  ``2012_HF/ref/<DOY>smoothed1.csv`` and ``tra/<DOY>smoothed1_tra.csv``; empty
  (all-NaN or all-zero) columns are leaves not sampled that day.
* Spectra columns are in the order of the traits workbook's sample columns. For
  Harvard Forest this is also the header order given in ``2012HFReadMe.txt``
  (checked here). spectra_db (and curated-leafspec before it) attached columns
  to samples in *alphabetical* order instead, which mislabeled every leaf.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

from leafrtm.db.adapters.core import DatasetTables, to_tables
from leafrtm.db.sources import yang2016
from leafrtm.db.sources.common import raw_dir

SITES = {
    "MV": {"year": 2011, "traits": "2011_MV_leaftraits.xlsx",
           "spectra": {"R": "2011-MV/*smoothed.csv"},
           # A/B/C = tree; U/L = upper/lower canopy; N/S = north/south side; leaf 1-3
           "sample": r"^(?P<tree>[ABC])(?P<canopy>[UL])(?P<aspect>[NS]?)(?P<leaf>\d)$"},
    "HF": {"year": 2012, "traits": "2012_HF_leaftraits.xlsx",
           "spectra": {"R": "2012_HF/ref/*smoothed1.csv", "T": "2012_HF/tra/*smoothed1_tra.csv"},
           # RO/RM/YB = species; tree number; U/L = upper/lower canopy; leaf number
           "sample": r"^(?P<species_code>RO|RM|YB)(?P<treeno>\d)(?P<canopy>[UL])(?P<leaf>\d)$"},
}
SHEETS = {"chla": "Chla", "chlb": "Chlb", "totchl": "TotChl", "car": "Car", "lma": "LMA",
          "totc": "TotC", "totn": "TotN"}


def _is_num(x) -> bool:
    return isinstance(x, (int, float, np.integer, np.floating)) and not pd.isna(x)


def _norm_name(x) -> str:
    return re.sub(r"[\s\-]", "", str(x))


def read_traits(path) -> tuple[list[str], pd.DataFrame]:
    """Returns (sample names in column order, long table doy x sample x trait columns)."""
    sheets = pd.read_excel(path, sheet_name=None, header=None)
    unknown = set(s.lower() for s in sheets) - set(SHEETS)
    if unknown:
        raise ValueError(f"{path.name}: unexpected sheets {sorted(unknown)}")
    ref = sheets[next(s for s in sheets if s.lower() == "chla")]
    names = [_norm_name(v) for v in ref.iloc[0, 1:] if pd.notna(v)]
    n = len(names)
    long = []
    for sname, df in sheets.items():
        trait = SHEETS[sname.lower()]
        is_doy = df.iloc[:, 0].map(_is_num)
        head = df[~is_doy].iloc[:1, 1:]  # header row (sample names) if present
        if len(head) and head.iloc[0].notna().any():
            h = [_norm_name(v) for v in head.iloc[0] if pd.notna(v)]
            if len(h) != n or not all(r.startswith(x) for r, x in zip(names, h)):
                raise ValueError(f"{path.name}/{sname}: header {h} does not match {names}")
        vals = df[is_doy]
        extra = vals.iloc[:, n + 1:]
        if extra.notna().any().any():
            raise ValueError(f"{path.name}/{sname}: values beyond the {n} sample columns")
        block = vals.iloc[:, 1:n + 1].apply(pd.to_numeric, errors="coerce")
        block.columns = names[:block.shape[1]]
        block.index = vals.iloc[:, 0].astype(int).rename("doy")
        long.append(block.stack(future_stack=True).rename(trait))
    traits = pd.concat(long, axis=1).rename_axis(["doy", "sample"]).reset_index()
    return names, traits


def read_spectra_file(path, names):
    d = pd.read_csv(path, header=None)
    d = d[d[0] != 0]  # trailing padding rows
    if d.shape[1] - 1 != len(names):
        raise ValueError(f"{path.name}: {d.shape[1] - 1} spectra columns for {len(names)} samples")
    wl = d[0].to_numpy(float)
    out = {}
    for j, s in enumerate(names, start=1):
        v = d[j].to_numpy(float)
        if np.isfinite(v).any() and np.nanmax(np.abs(v)) > 0:  # empty column: no sample
            out[s] = (j, v)
    return wl, out


def load_site(code: str, root):
    spec = SITES[code]
    names, traits = read_traits(root / spec["traits"])
    doys = sorted(traits["doy"].unique())
    if code == "HF":  # the readme's header order must be the workbook's sample order
        hdr = [_norm_name(x) for x in (root / "2012HFReadMe.txt").read_text().split("\n")[0].split()]
        if hdr != names:
            raise ValueError(f"2012HFReadMe.txt headers {hdr} differ from the workbook {names}")
    rows, spectra = [], []
    for q, pattern in spec["spectra"].items():
        files = sorted(root.glob(pattern))
        if code == "MV":  # NNsmoothed.csv = NN-th sampling day
            idx = [int(re.match(r"(\d+)smoothed\.csv$", f.name).group(1)) for f in files]
            if sorted(idx) != list(range(1, len(doys) + 1)):
                raise ValueError(f"MV spectra files {idx} do not match {len(doys)} sampling days")
            file_doy = {f: doys[i - 1] for f, i in zip(files, idx)}
        else:
            file_doy = {f: int(re.match(r"(\d+)smoothed1", f.name).group(1)) for f in files}
            if set(file_doy.values()) != set(doys):
                raise ValueError(f"HF {q} files {sorted(file_doy.values())} != trait days {doys}")
        for f, doy in file_doy.items():
            wl, cols = read_spectra_file(f, names)
            for s, (j, v) in cols.items():
                rows.append({"sample": s, "doy": doy, "quantity": q, "file": f.name, "file_column": j})
                spectra.append(pd.Series(v, index=wl))
    meta = pd.DataFrame(rows)
    parts = meta["sample"].str.extract(spec["sample"])
    if parts.isna().all(axis=1).any():
        raise ValueError(f"{code}: unparsed sample names {sorted(meta['sample'][parts.isna().all(axis=1)])}")
    meta = pd.concat([meta, parts], axis=1)
    if code == "MV":
        meta["species_code"] = "RO"  # all trees are red oak (workbook note)
        meta["aspect"] = meta["aspect"].replace("", None)
    else:
        meta["tree"] = meta["species_code"] + meta["treeno"]
        meta = meta.drop(columns="treeno")
    meta["site"], meta["year"] = code, spec["year"]
    meta["year_doy"] = meta["year"].astype(str) + "-" + meta["doy"].astype(str).str.zfill(3)
    # same IDs as the spectra_db version (referenced by earlier fit sets)
    meta["observation_id"] = ("yang_pheno|" + meta["sample"] + "|" + meta["doy"].astype(str) + "|"
                              + meta["year"].astype(str))
    with_spec = set(zip(meta["doy"], meta["sample"]))
    tcols = [c for c in traits.columns if c not in ("doy", "sample")]
    has_trait = traits[tcols].notna().any(axis=1)
    orphan = [(d, s) for d, s in zip(traits["doy"][has_trait], traits["sample"][has_trait])
              if (d, s) not in with_spec]
    meta = meta.merge(traits, on=["doy", "sample"], how="left")
    return meta, spectra, orphan


def load(cfg: dict) -> DatasetTables:
    root = raw_dir("yang2016")
    if not all((root / s["traits"]).exists() for s in SITES.values()) or \
            not (root / "2012HFReadMe.txt").exists():
        yang2016.fetch_all()
    metas, spectra, orphans = [], [], {}
    for code in SITES:
        m, s, o = load_site(code, root)
        metas.append(m)
        spectra += s
        orphans[code] = o
    meta = pd.concat(metas, ignore_index=True)
    wl = np.unique(np.concatenate([s.index.to_numpy() for s in spectra]))
    vals = np.full((len(spectra), len(wl)), np.nan)
    for i, s in enumerate(spectra):
        vals[i, np.searchsorted(wl, s.index.to_numpy())] = s.to_numpy()
    info = {"title": "Leaf traits and spectral data from Yang et al., 2016 RSE",
            "url": f"https://github.com/{yang2016.REPO}/tree/{yang2016.COMMIT}",
            "doi": cfg.get("doi"), "citation": cfg.get("citation"), "license": None}
    out = to_tables(cfg, meta, wl, vals, info)
    # leaf-days with trait values but no spectrum (not in the output)
    out.report["trait_samples_without_spectra"] = {k: [f"{s}@{d}" for d, s in v] for k, v in orphans.items()}
    return out
