"""Config-driven conversion of a wide "one row per spectrum" table into the
normalized schema (``leafrtm.db.schema``).

Every source module produces ``(meta, wavelength, values)``: a metadata
DataFrame with one row per spectrum, the wavelength vector, and an
(n_spectra x n_bands) array. A per-dataset YAML config (``configs/*.yaml``;
format documented in ``configs/README.md``) says how metadata columns map to
IDs, species, dates, sites, treatments, enums and traits.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from functools import cache
from importlib.resources import files
from pathlib import Path

import numpy as np
import pandas as pd
import pint
import yaml

from leafrtm.db import schema

UREG = pint.UnitRegistry()
UREG.define("permille = 0.001 = ‰")
CONFIG_DIR = Path(__file__).parent / "configs"


@cache
def trait_vocab() -> dict:
    with (files("leafrtm.db") / "vocab" / "traits.yaml").open() as f:
        return yaml.safe_load(f)


@cache
def load_config(dataset_id: str) -> dict:
    with open(CONFIG_DIR / f"{dataset_id}.yaml") as f:
        cfg = yaml.safe_load(f)
    if cfg.get("dataset_id") != dataset_id:
        raise ValueError(f"{dataset_id}.yaml declares dataset_id={cfg.get('dataset_id')}")
    return cfg


def list_configs(source: str | None = None, include_only=True) -> list[str]:
    out = []
    for p in sorted(CONFIG_DIR.glob("*.yaml")):
        cfg = load_config(p.stem)
        if include_only and not cfg.get("include", True):
            continue
        if source and cfg.get("source") != source:
            continue
        out.append(p.stem)
    return out


@dataclass
class DatasetTables:
    """Normalized tables for one dataset (before cross-dataset merging)."""

    tables: dict[str, pd.DataFrame] = field(default_factory=dict)
    grids: dict[str, tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)
    raw_species: pd.DataFrame | None = None  # entity_id -> raw species fields
    report: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Column mapping helpers
# ---------------------------------------------------------------------------


def _norm(v):
    if isinstance(v, str):
        return v.strip()
    return v


def map_column(meta: pd.DataFrame, spec, default=None, enum=None) -> pd.Series:
    """Resolve a config entry that is a constant or {column, map, default}."""
    if spec is None:
        return pd.Series(default, index=meta.index, dtype=object)
    if not isinstance(spec, dict):
        return pd.Series(spec, index=meta.index, dtype=object)
    if "columns" in spec:  # concatenate several columns
        s = meta[spec["columns"]].astype(str).agg(spec.get("sep", " ").join, axis=1)
    else:
        s = meta[spec["column"]]
    s = s.map(_norm)
    if "map" in spec:
        m = {(_norm(k) if not isinstance(k, str) else k.strip()): v for k, v in spec["map"].items()}
        mapped = s.map(lambda x: m.get(x, m.get(str(x), m.get(str(x).lower(), np.nan))))
        if "default" in spec:
            mapped = mapped.where(mapped.notna(), spec["default"])
        s = mapped
    elif "regex" in spec:  # list of [pattern, value]
        out = pd.Series(spec.get("default", np.nan), index=meta.index, dtype=object)
        for pat, val in reversed(spec["regex"]):
            hit = s.astype(str).str.contains(pat, case=False, regex=True, na=False)
            out[hit] = val
        s = out
    s = s.where(pd.notna(s), None)
    if enum is not None:
        bad = set(s.dropna()) - set(enum)
        if bad:
            raise ValueError(f"Values {sorted(map(str, bad))} not in enum {enum}")
    return s


def _id_str(x) -> str:
    """Stable string form of an ID value (1.0 -> "1"; missing -> "NA")."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "NA"
    if isinstance(x, float) and x.is_integer():
        return str(int(x))
    return str(x).strip()


def _id_cols(meta, cols):
    if not cols:
        return pd.Series("", index=meta.index)
    if isinstance(cols, str):
        cols = [cols]
    return meta[cols].map(_id_str).agg("|".join, axis=1)


def parse_dates(meta, spec) -> tuple[pd.Series, pd.Series]:
    """Returns (date, year)."""
    idx = meta.index
    if spec is None:
        return pd.Series(pd.NaT, index=idx), pd.Series(pd.NA, index=idx, dtype="Int64")
    if not isinstance(spec, dict):  # constant year or date
        s = str(spec)
        if re.fullmatch(r"\d{4}", s):
            return pd.Series(pd.NaT, index=idx), pd.Series(int(s), index=idx, dtype="Int64")
        d = pd.Timestamp(s)
        return pd.Series(d, index=idx), pd.Series(d.year, index=idx, dtype="Int64")
    if "year" in spec and "column" not in spec:
        y = pd.to_numeric(meta[spec["year"]], errors="coerce").astype("Int64")
        return pd.Series(pd.NaT, index=idx), y
    raw = meta[spec["column"]]
    if spec.get("format") == "excel":
        d = pd.to_datetime(pd.to_numeric(raw, errors="coerce"), unit="D", origin="1899-12-30")
    elif "format" in spec:
        d = pd.to_datetime(raw.astype(str).str.strip(), format=spec["format"], errors="coerce")
    else:
        d = pd.to_datetime(raw, errors="coerce")
    if "year" in spec:
        y = pd.to_numeric(meta[spec["year"]], errors="coerce").astype("Int64")
        y = y.where(y.notna(), d.dt.year.astype("Int64"))
    else:
        y = d.dt.year.astype("Int64")
    return d, y


def convert(values: pd.Series, unit: str | None, trait: str, factor: float = 1.0) -> pd.Series:
    target = trait_vocab()[trait]["unit"]
    v = pd.to_numeric(values, errors="coerce") * factor
    if unit is None or unit == target:
        return v
    f = UREG.Quantity(1.0, unit).to(target).magnitude
    return v * f


# ---------------------------------------------------------------------------
# Spectra
# ---------------------------------------------------------------------------


def grid_id(wavelength, segment_start) -> str:
    h = hashlib.sha1(np.asarray(wavelength, np.float64).tobytes()
                     + np.asarray(segment_start, bool).tobytes()).hexdigest()[:16]
    return f"g{h}"


def segment_starts(wl, joins):
    s = np.zeros(len(wl), bool)
    if len(wl):
        s[0] = True
    for j in joins or []:
        i = np.searchsorted(wl, j)
        if 0 < i < len(wl):
            s[i] = True
    return s


def scale_spectra(values: np.ndarray, quantity: np.ndarray, mode: str) -> tuple[np.ndarray, np.ndarray]:
    """Convert percent to fraction. Returns (values, was_percent)."""
    v = values.copy()
    if mode == "fraction":
        return v, np.zeros(len(v), bool)
    med = np.nanmedian(np.where(np.isfinite(v), v, np.nan), axis=1)
    if mode == "percent":
        pct = np.ones(len(v), bool)
    else:  # auto, per spectrum: R/T/A medians in fraction are < 1
        pct = med > 1.5
    pct &= np.isin(quantity, ["R", "T", "A"])
    v[pct] /= 100.0
    return v, pct


# ---------------------------------------------------------------------------
# Main conversion
# ---------------------------------------------------------------------------


def to_tables(cfg: dict, meta: pd.DataFrame, wavelength, values, dataset_meta: dict) -> DatasetTables:
    ds = cfg["dataset_id"]
    out = DatasetTables()
    meta = meta.reset_index(drop=True).copy()
    values = np.array(values, dtype=np.float64)  # copy: modified in place
    wavelength = np.asarray(wavelength, dtype=np.float64)
    rep = {"rows_in": len(meta)}

    sentinels = cfg.get("sentinels", [-999, -9999, -99999])
    meta = meta.replace(sentinels, np.nan)
    values[np.isin(values, sentinels)] = np.nan

    for name, expr in cfg.get("columns", {}).items():  # derived metadata columns
        meta[name] = meta.eval(expr, engine="python")

    for q in np.atleast_1d(cfg.get("filter", [])):
        keep = meta.eval(q, engine="python").to_numpy(bool)
        meta, values = meta[keep].reset_index(drop=True), values[keep]
    rep["rows_after_filter"] = len(meta)

    sp = cfg.get("spectra", {})
    quantity = map_column(meta, sp.get("quantity", "R"), enum=schema.QUANTITY + [None])
    keep = quantity.notna().to_numpy()
    meta, values, quantity = meta[keep].reset_index(drop=True), values[keep], quantity[keep].reset_index(drop=True)
    rep["rows_with_quantity"] = len(meta)

    lo, hi = sp.get("wavelength_range", [-np.inf, np.inf])
    inr = (wavelength >= lo) & (wavelength <= hi)
    wavelength, values = wavelength[inr], values[:, inr]
    order = np.argsort(wavelength, kind="stable")
    wavelength, values = wavelength[order], values[:, order]
    if np.any(np.diff(wavelength) == 0):  # duplicated bands: average
        uniq, inv = np.unique(wavelength, return_inverse=True)
        acc = np.zeros((len(values), len(uniq)))
        cnt = np.zeros((len(values), len(uniq)))
        fin = np.isfinite(values)
        np.add.at(acc.T, inv, np.where(fin, values, 0).T)
        np.add.at(cnt.T, inv, fin.T.astype(float))
        with np.errstate(invalid="ignore"):
            values = acc / cnt
        wavelength = uniq
    values, pct = scale_spectra(values, quantity.to_numpy(), sp.get("scale", "auto"))
    rep["percent_scaled"] = int(pct.sum())

    # --- IDs -----------------------------------------------------------------
    ids = cfg["ids"]
    ent_local = _id_cols(meta, ids["entity"])
    time_local = _id_cols(meta, ids.get("time"))
    date, year = parse_dates(meta, cfg.get("date"))
    entity_id = ds + ":" + ent_local
    obs_local = ent_local + np.where(time_local != "", "@" + time_local, "")
    observation_id = ds + ":" + obs_local

    # --- datasets / sites / plots ---------------------------------------------
    out.tables["datasets"] = pd.DataFrame([{
        "dataset_id": ds, "title": cfg.get("title") or dataset_meta.get("title") or ds,
        "source": cfg["source"], "source_id": cfg.get("source_id", ds),
        "url": cfg.get("url") or dataset_meta.get("url"),
        "doi": cfg.get("doi") or dataset_meta.get("doi"),
        "citation": cfg.get("citation") or dataset_meta.get("citation"),
        "license": cfg.get("license") or dataset_meta.get("license"),
    }])
    site_spec = cfg.get("site", {"name": ds})
    if isinstance(site_spec, dict) and "column" in site_spec:
        site_name = map_column(meta, site_spec).astype(object).where(lambda s: s.notna(), "unknown")
    else:
        site_name = pd.Series(site_spec.get("name", ds) if isinstance(site_spec, dict) else site_spec,
                              index=meta.index)
    site_id = ds + ":" + site_name.astype(str)
    lat = pd.to_numeric(meta[cfg["lat"]], errors="coerce") if isinstance(cfg.get("lat"), str) else \
        pd.Series(cfg.get("lat", np.nan), index=meta.index, dtype=float)
    lon = pd.to_numeric(meta[cfg["lon"]], errors="coerce") if isinstance(cfg.get("lon"), str) else \
        pd.Series(cfg.get("lon", np.nan), index=meta.index, dtype=float)
    if cfg.get("crs"):
        from pyproj import Transformer  # optional; only for projected coordinates
        tr = Transformer.from_crs(cfg["crs"], "EPSG:4326", always_xy=True)
        lon_, lat_ = tr.transform(lon.to_numpy(), lat.to_numpy())
        lat, lon = pd.Series(lat_, index=meta.index), pd.Series(lon_, index=meta.index)
    bad = (lat.abs() > 90) | (lon.abs() > 180)
    lat[bad], lon[bad] = np.nan, np.nan
    plot_spec = cfg.get("plot")
    plot_name = map_column(meta, plot_spec) if plot_spec else pd.Series(None, index=meta.index, dtype=object)
    plot_id = pd.Series(np.where(plot_name.notna(), site_id + ":" + plot_name.astype(str), None),
                        index=meta.index, dtype=object)
    geo = pd.DataFrame({"site_id": site_id, "site": site_name, "plot_id": plot_id,
                        "plot": plot_name, "lat": lat, "lon": lon})
    sites = geo.groupby("site_id", as_index=False).agg(name=("site", "first"), lat=("lat", "mean"),
                                                       lon=("lon", "mean"))
    sites["dataset_id"] = ds
    sites["ecosystem"] = cfg.get("ecosystem")
    out.tables["sites"] = sites[["site_id", "dataset_id", "name", "lat", "lon", "ecosystem"]]
    g = geo.dropna(subset=["plot_id"])
    out.tables["plots"] = g.groupby("plot_id", as_index=False).agg(
        site_id=("site_id", "first"), name=("plot", "first"), lat=("lat", "mean"), lon=("lon", "mean"))

    # --- entities --------------------------------------------------------------
    sp_cfg = cfg.get("species", {})
    raw = pd.DataFrame({"entity_id": entity_id})
    if "constant" in sp_cfg:
        raw["name"] = sp_cfg["constant"]
    elif "name" in sp_cfg:
        cols = sp_cfg["name"] if isinstance(sp_cfg["name"], list) else [sp_cfg["name"]]
        parts = meta[cols].map(lambda x: x.strip() if isinstance(x, str) else x)
        raw["name"] = parts.apply(lambda r: " ".join(str(v) for v in r if pd.notna(v) and str(v)),
                                  axis=1).replace("", None)
        if "map" in sp_cfg:
            raw["name"] = raw["name"].map(lambda x: sp_cfg["map"].get(x, x))
    if "code" in sp_cfg:
        raw["code"] = meta[sp_cfg["code"]].map(lambda x: str(x).strip() if pd.notna(x) else None)
        raw["code_type"] = sp_cfg.get("code_type", "usda")
        if "code_map" in sp_cfg:
            raw["name"] = raw["code"].map(sp_cfg["code_map"]).where(
                raw["code"].isin(sp_cfg["code_map"]), raw.get("name"))
    out.raw_species = raw.drop_duplicates("entity_id")

    ent = pd.DataFrame({
        "entity_id": entity_id, "dataset_id": ds, "plot_id": plot_id,
        "species_id": None,  # resolved across datasets in harmonize.py
        "plant_id": (ds + ":" + _id_cols(meta, ids["plant"])) if ids.get("plant") else None,
        "canopy_position": map_column(meta, cfg.get("canopy_position"), enum=schema.CANOPY_POSITION + [None]),
        "sun_shade": map_column(meta, cfg.get("sun_shade"), enum=schema.SUN_SHADE + [None]),
        "leaf_age": map_column(meta, cfg.get("leaf_age"), enum=schema.LEAF_AGE + [None]),
    })
    out.tables["entities"] = ent.groupby("entity_id", as_index=False).first()

    # --- observations ------------------------------------------------------------
    obs = pd.DataFrame({
        "observation_id": observation_id, "entity_id": entity_id, "dataset_id": ds,
        "date": date, "year": year,
        "sample_state": map_column(meta, sp.get("sample_state", "fresh"), enum=schema.SAMPLE_STATE + [None]),
        "target_config": map_column(meta, sp.get("target_config", "single_leaf"),
                                    enum=schema.TARGET_CONFIG + [None]),
    })
    if obs["sample_state"].isna().any():
        raise ValueError(f"{ds}: unmapped sample_state values")
    obs["target_config"] = obs["target_config"].fillna("unknown")
    dup = obs.groupby("observation_id")[["sample_state"]].nunique()
    if (dup["sample_state"] > 1).any():
        raise ValueError(f"{ds}: observations mix sample states; add the state column to ids.time")
    out.tables["observations"] = obs.groupby("observation_id", as_index=False).first()

    # --- treatments --------------------------------------------------------------
    trt_rows, link_rows = [], []
    for t in cfg.get("treatments", []):
        lvl = map_column(meta, {k: v for k, v in t.items() if k in ("column", "columns", "map", "sep", "regex")})
        controls = {str(c) for c in np.atleast_1d(t.get("control", []))}
        for level in sorted(set(lvl.dropna().astype(str))):
            tid = f"{ds}:{t['factor']}={level}"
            trt_rows.append({"treatment_id": tid, "dataset_id": ds, "factor": t["factor"],
                             "level": level,
                             "is_control": (level in controls) if controls else pd.NA})
        has = lvl.notna()
        link_rows.append(pd.DataFrame({
            "observation_id": observation_id[has],
            "treatment_id": ds + ":" + t["factor"] + "=" + lvl[has].astype(str)}))
    out.tables["treatments"] = pd.DataFrame(
        trt_rows, columns=["treatment_id", "dataset_id", "factor", "level", "is_control"]).astype(
        {"is_control": "boolean"})
    out.tables["observation_treatments"] = (pd.concat(link_rows).drop_duplicates() if link_rows else
                                            pd.DataFrame(columns=["observation_id", "treatment_id"]))

    # --- instruments ---------------------------------------------------------------
    instrument = map_column(meta, sp.get("instrument"), default=None)
    foreoptic = map_column(meta, sp.get("foreoptic", "unknown"), enum=schema.FOREOPTIC + [None])
    foreoptic = foreoptic.fillna("unknown")

    # --- spectra -------------------------------------------------------------------
    joins = sp.get("segments", [])
    rows = []
    rep_idx = pd.Series(observation_id).groupby([observation_id, quantity]).cumcount()
    finite = np.isfinite(values)
    keys = [m.tobytes() for m in finite]
    by_mask = {}
    for i, k in enumerate(keys):
        by_mask.setdefault(k, []).append(i)
    for k, idx in by_mask.items():
        m = finite[idx[0]]
        if m.sum() < 10:
            continue
        wl = wavelength[m]
        seg = segment_starts(wl, joins)
        gid = grid_id(wl, seg)
        out.grids[gid] = (wl, seg)
        for i in idx:
            rows.append({
                "spectra_id": f"{observation_id[i]}#{quantity[i]}{rep_idx[i] + 1}",
                "observation_id": observation_id[i], "quantity": quantity[i],
                "instrument_id": instrument[i], "foreoptic": foreoptic[i], "grid_id": gid,
                "qc_flags": None, "values": values[i, m].astype(np.float32)})
    spectra = pd.DataFrame(rows)
    if spectra["spectra_id"].duplicated().any():
        raise ValueError(f"{ds}: duplicated spectra ids")
    out.tables["spectra"] = spectra
    rep["spectra"] = len(spectra)

    # --- traits ----------------------------------------------------------------------
    trows = []
    for t in cfg.get("traits", []):
        if t["trait"] not in trait_vocab():
            raise ValueError(f"{ds}: unknown trait {t['trait']}")
        v = convert(meta[t["column"]], t.get("unit"), t["trait"], t.get("factor", 1.0))
        level = t.get("level", "leaf")
        df = pd.DataFrame({
            "observation_id": observation_id if level == "leaf" else None,
            "entity_id": entity_id, "trait": t["trait"], "value": v,
            "unit": trait_vocab()[t["trait"]]["unit"], "trait_level": level,
            "original_name": t["column"], "original_value": pd.to_numeric(meta[t["column"]], errors="coerce"),
            "original_unit": t.get("unit"), "method": t.get("method"), "derived": False})
        df = df[np.isfinite(df["value"])]
        # one value per observation (replicate rows repeat trait values)
        key = "observation_id" if level == "leaf" else "entity_id"
        trows.append(df.groupby([key, "trait"], as_index=False, dropna=False).first())
    out.tables["traits"] = (pd.concat(trows, ignore_index=True) if trows else
                            pd.DataFrame(columns=list(schema.Traits.to_schema().columns)))
    rep["traits"] = len(out.tables["traits"])
    rep["observations"] = len(out.tables["observations"])
    out.report = rep
    return out
