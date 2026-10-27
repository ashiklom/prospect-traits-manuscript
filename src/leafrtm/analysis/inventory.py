"""Summaries of the harmonized database for the website's Data pages.

    python -m leafrtm.analysis.inventory --out results/pilot/site_data
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from leafrtm.db.harmonize import load, processed_dir


def inventory(db: dict) -> dict[str, pd.DataFrame]:
    obs, sp, tr = db["observations"], db["spectra"], db["traits"]
    q = sp.groupby("observation_id")["quantity"].agg(lambda s: "+".join(sorted(set(s))))
    o = obs.join(q.rename("quantities"), on="observation_id").merge(
        db["entities"][["entity_id", "species_id"]], on="entity_id", how="left")
    fo = sp.groupby("observation_id")["foreoptic"].first()
    o = o.join(fo.rename("foreoptic"), on="observation_id")
    sites = db["sites"].merge(o.groupby("dataset_id").size().rename("n_obs"), on="dataset_id")
    ds = db["datasets"].set_index("dataset_id")
    g = o.groupby("dataset_id")
    ds["n_observations"] = g.size()
    ds["n_spectra"] = sp.assign(d=sp["observation_id"].str.split(":").str[0]).groupby("d").size()
    ds["n_species"] = g["species_id"].nunique()
    ds["quantities"] = g["quantities"].agg(lambda s: ", ".join(sorted(set(s.dropna()))))
    ds["sample_states"] = g["sample_state"].agg(lambda s: ", ".join(sorted(set(s))))
    ds["foreoptic"] = g["foreoptic"].agg(lambda s: ", ".join(sorted(set(s.dropna()))))
    trd = tr.assign(d=tr["entity_id"].str.split(":").str[0])
    ds["traits"] = trd[~trd["derived"]].groupby("d")["trait"].agg(lambda s: ", ".join(sorted(set(s))))
    trt = db["treatments"].groupby("dataset_id")["factor"].agg(lambda s: ", ".join(sorted(set(s))))
    ds["treatments"] = trt
    ds["lat"] = db["sites"].groupby("dataset_id")["lat"].mean()
    ds["lon"] = db["sites"].groupby("dataset_id")["lon"].mean()
    coverage = trd.groupby(["d", "trait"]).size().rename("n").reset_index().rename(
        columns={"d": "dataset_id"})
    states = o.groupby(["sample_state", "quantities", "foreoptic"], dropna=False).size().rename(
        "n").reset_index()
    flags = sp.assign(d=sp["observation_id"].str.split(":").str[0],
                      flag=sp["qc_flags"].fillna("ok").str.split(";")).explode("flag")
    qc = flags.groupby(["d", "flag"]).size().rename("n").reset_index().rename(columns={"d": "dataset_id"})
    species = db["species"].merge(o.groupby("species_id").agg(
        n_observations=("observation_id", "size"),
        datasets=("dataset_id", lambda s: ", ".join(sorted(set(s))))).reset_index(), on="species_id")
    return {"inventory_datasets": ds.reset_index(), "inventory_sites": sites,
            "inventory_trait_coverage": coverage, "inventory_states": states,
            "inventory_qc": qc, "inventory_species": species,
            "inventory_trait_summary": tr.groupby(["trait", "unit"])["value"].describe().reset_index()}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for name, df in inventory(load()).items():
        df.to_parquet(out / f"{name}.parquet", index=False)
    rep = json.loads((processed_dir() / "report.json").read_text())
    (out / "harmonize_report.json").write_text(json.dumps(rep, indent=1, default=str))
    review = processed_dir() / "species_review.csv"
    if review.exists():
        pd.read_csv(review).to_parquet(out / "species_review.parquet", index=False)


if __name__ == "__main__":
    main()
