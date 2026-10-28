"""Merge per-dataset tables into the harmonized database.

    python -m leafrtm.db.harmonize [--rebuild] [dataset_id ...]

Steps: convert each included dataset (cached under ``data/interim``),
concatenate, resolve species, derive traits, screen implausible trait
values, flag spectra (QC), validate the schema and foreign keys, and write
``data/processed/db/<table>.parquet`` plus ``report.json`` and
``species_review.csv`` (names needing manual review).
"""

from __future__ import annotations

import argparse
import json
from importlib.resources import files

import numpy as np
import pandas as pd
import yaml

from leafrtm.db import qc, schema
from leafrtm.db import species as sp
from leafrtm.db.adapters.core import list_configs, trait_vocab
from leafrtm.db.build import convert, interim_dir, write
from leafrtm.db.sources.common import data_dir

TABLE_NAMES = ["datasets", "sites", "plots", "treatments", "entities", "observations",
               "observation_treatments", "spectra", "traits"]


def processed_dir():
    d = data_dir() / "processed" / "db"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _load_interim(ds, rebuild):
    d = interim_dir(ds)
    if rebuild or not (d / "spectra.parquet").exists():
        write(ds, convert(ds))
    out = {t: pd.read_parquet(d / f"{t}.parquet") for t in TABLE_NAMES}
    out["wavelength_grids"] = pd.read_parquet(d / "wavelength_grids.parquet")
    rs = pd.read_parquet(d / "raw_species.parquet")
    rs["dataset_id"] = ds
    out["raw_species"] = rs
    return out


# ---------------------------------------------------------------------------
# Derived traits
# ---------------------------------------------------------------------------


def derive_traits(traits: pd.DataFrame) -> pd.DataFrame:
    """Add traits computable from others measured on the same observation.

    EWT = LMA * W / (1 - W) with W the water fraction of fresh mass;
    LWC from EWT and LMA; N per area <-> per mass via LMA; Chl a+b = a + b.
    """
    leaf = traits[traits["trait_level"] == "leaf"]
    wide = leaf.pivot_table(index=["observation_id", "entity_id"], columns="trait",
                            values="value", aggfunc="mean")
    new = {}

    def has(*cols):
        return all(c in wide for c in cols)

    if has("lma", "lwc"):
        w = wide["lwc"] / 100
        new["ewt"] = wide["lma"] * w / (1 - w)
    if has("lma", "ewt"):
        new["lwc"] = 100 * wide["ewt"] / (wide["ewt"] + wide["lma"])
    if has("lma", "n_mass"):
        new["n_area"] = wide["n_mass"] / 100 * wide["lma"]
    if has("lma", "n_area"):
        new["n_mass"] = 100 * wide["n_area"] / wide["lma"]
    if has("chla_area", "chlb_area"):
        new["chl_area"] = wide["chla_area"] + wide["chlb_area"]
    if has("lma", "c_mass"):
        new["c_area"] = wide["c_mass"] / 100 * wide["lma"]
    rows = []
    for trait, v in new.items():
        if trait in wide:
            v = v[wide[trait].isna()]  # only where not measured
        v = v.dropna()
        if v.empty:
            continue
        df = v.rename("value").reset_index()
        df["trait"] = trait
        rows.append(df)
    if not rows:
        return traits
    d = pd.concat(rows, ignore_index=True)
    d = d[np.isfinite(d["value"].astype(float))]
    d["unit"] = d["trait"].map(lambda t: trait_vocab()[t]["unit"])
    d["trait_level"] = "leaf"
    d["original_name"] = "derived"
    d["original_value"] = np.nan
    d["original_unit"] = None
    d["method"] = "derived"
    d["derived"] = True
    return pd.concat([traits, d[traits.columns]], ignore_index=True)


# ---------------------------------------------------------------------------
# Species and instruments
# ---------------------------------------------------------------------------


def species_tables(raw: pd.DataFrame):
    res = sp.resolve(raw)
    res["pass1_name"] = res["resolved"]
    res["resolved"] = sp.accepted_names(res["resolved"])
    res["species_id"] = res["resolved"].map(sp.species_id)
    names = res.dropna(subset=["resolved"]).drop_duplicates("species_id")
    species = pd.DataFrame({"species_id": names["species_id"], "scientific_name": names["resolved"]})
    species["genus"] = species["scientific_name"].str.split(" ").str[0]
    gb = sp._read_vocab(sp.GBIF)
    if len(gb):
        # GBIF matched the pass-1 names; carry each family to its accepted name
        fam = gb.drop_duplicates("name").set_index("name")["family"]
        by_acc = (res.dropna(subset=["resolved"]).assign(family=lambda d: d["pass1_name"].map(fam))
                  .dropna(subset=["family"]).drop_duplicates("resolved").set_index("resolved")["family"])
        species["family"] = species["scientific_name"].map(by_acc)
    else:
        u = sp.usda().drop_duplicates("name").set_index("name")
        species["family"] = species["scientific_name"].map(u["family"])
    attrs = sp._read_vocab(sp.ATTRIBUTES, dtype=str)
    for col in ("growth_form", "leaf_type", "phenology", "ps_pathway"):
        species[col] = None
        if len(attrs) and col in attrs:
            a = attrs.dropna(subset=[col])
            by_sp = dict(zip(a["taxon"], a[col]))
            species[col] = species["scientific_name"].map(by_sp)
            by_genus = {k: v for k, v in by_sp.items() if " " not in k}
            species[col] = species[col].fillna(species["genus"].map(by_genus))
    return species.reset_index(drop=True), res[["entity_id", "species_id"]], res


def instruments_table(ids) -> pd.DataFrame:
    with (files("leafrtm.db") / "vocab" / "instruments.yaml").open() as f:
        voc = yaml.safe_load(f)
    rows = []
    for i in sorted(set(ids) - {None}):
        v = voc.get(i, {})
        rows.append({"instrument_id": i, "manufacturer": v.get("manufacturer"),
                     "model": v.get("model"), "fwhm_nm": v.get("fwhm_nm")})
    return pd.DataFrame(rows, columns=["instrument_id", "manufacturer", "model", "fwhm_nm"])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def harmonize(dataset_ids=None, rebuild=False) -> dict:
    ids = dataset_ids or list_configs()
    parts = [_load_interim(ds, rebuild) for ds in ids]
    db = {t: pd.concat([p[t] for p in parts], ignore_index=True) for t in TABLE_NAMES}
    db["wavelength_grids"] = pd.concat([p["wavelength_grids"] for p in parts]).drop_duplicates(
        "grid_id").reset_index(drop=True)
    raw = pd.concat([p["raw_species"] for p in parts], ignore_index=True)
    db["datasets"]["source_id"] = db["datasets"]["source_id"].map(
        lambda v: "; ".join(map(str, v)) if not isinstance(v, str) else v)
    for t in ("sites", "plots"):  # numeric plot/site labels (e.g. plot 117) are names, not numbers
        db[t]["name"] = db[t]["name"].map(lambda v: v if pd.isna(v) else str(v)).astype(object)
    # datasets differ in date types (tz-aware, naive, missing): normalize to naive UTC dates
    db["observations"]["date"] = pd.to_datetime(db["observations"]["date"], errors="coerce",
                                                utc=True).dt.tz_localize(None)
    for c in ("value", "original_value"):  # empty per-dataset tables concat as object
        db["traits"][c] = pd.to_numeric(db["traits"][c], errors="coerce")
    db["traits"]["derived"] = db["traits"]["derived"].astype(bool)
    report = {"datasets": len(ids)}

    species, ent_sp, res = species_tables(raw)
    db["species"] = species
    db["entities"] = db["entities"].drop(columns="species_id").merge(ent_sp, on="entity_id", how="left")
    report["species"] = {"n": len(species), "resolution": res["resolution"].value_counts().to_dict()}
    sp.review_table(res).to_csv(processed_dir() / "species_review.csv", index=False)

    tr = derive_traits(db["traits"])
    tr, dropped = qc.screen_traits(tr)
    db["traits"] = tr
    report["traits_dropped_implausible"] = dropped.groupby(["trait"]).size().to_dict()
    report["traits_derived"] = tr[tr["derived"]].groupby("trait").size().to_dict()

    grids = {g: np.asarray(w) for g, w in zip(db["wavelength_grids"]["grid_id"],
                                             db["wavelength_grids"]["wavelength"])}
    s = db["spectra"]
    flags = [qc.spectral_flags(np.asarray(v, float), grids[g], q)
             for v, g, q in zip(s["values"], s["grid_id"], s["quantity"])]
    rep_out = qc.replicate_outliers(s, grids)
    s["qc_flags"] = [";".join(f + (["replicate_outlier"] if o else [])) or None
                     for f, o in zip(flags, rep_out)]
    report["qc_flags"] = s["qc_flags"].str.split(";").explode().value_counts().to_dict()
    db["instruments"] = instruments_table(s["instrument_id"].dropna().unique())

    for t in ("sites", "plots", "treatments", "observation_treatments"):
        key = {"sites": "site_id", "plots": "plot_id", "treatments": "treatment_id"}.get(t)
        db[t] = db[t].drop_duplicates(key) if key else db[t].drop_duplicates()
    schema.validate(db)
    report["counts"] = {t: len(v) for t, v in db.items()}
    report["by_dataset"] = (db["observations"].merge(
        s.groupby("observation_id")["quantity"].agg(lambda q: "".join(sorted(set(q)))).rename("q"),
        on="observation_id").groupby(["dataset_id", "sample_state", "q"]).size()
        .reset_index(name="n").to_dict("records"))
    out = processed_dir()
    for t, df in db.items():
        df.to_parquet(out / f"{t}.parquet", index=False)
    (out / "report.json").write_text(json.dumps(report, indent=1, default=str))
    return report


def load(tables=None) -> dict[str, pd.DataFrame]:
    d = processed_dir()
    names = tables or [p.stem for p in d.glob("*.parquet")]
    return {t: pd.read_parquet(d / f"{t}.parquet") for t in names}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebuild", action="store_true")
    ap.add_argument("dataset_ids", nargs="*")
    a = ap.parse_args(argv)
    rep = harmonize(a.dataset_ids or None, a.rebuild)
    print(json.dumps({k: v for k, v in rep.items() if k != "by_dataset"}, indent=1, default=str))


if __name__ == "__main__":
    main()
