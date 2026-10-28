"""Build normalized per-dataset tables from adapter configs.

    python -m leafrtm.db.build check <dataset_id> [...]   # convert + report, no output
    python -m leafrtm.db.build dataset <dataset_id> [...]  # write data/interim/<dataset_id>/
"""

from __future__ import annotations

import argparse
import importlib
import json
import sys

import numpy as np
import pandas as pd

from leafrtm.db import schema
from leafrtm.db.adapters.core import DatasetTables, load_config
from leafrtm.db.sources.common import data_dir


def interim_dir(dataset_id):
    d = data_dir() / "interim" / dataset_id
    d.mkdir(parents=True, exist_ok=True)
    return d


def convert(dataset_id: str) -> DatasetTables:
    cfg = load_config(dataset_id)
    if not cfg.get("include", True):
        raise ValueError(f"{dataset_id} is excluded: {cfg.get('exclude_reason')}")
    loader = importlib.import_module(f"leafrtm.db.adapters.{cfg['source']}")
    return loader.load(cfg)


def summarize(dt: DatasetTables) -> dict:
    t = dt.tables
    sp = t["spectra"]
    rep = dict(dt.report)
    rep["quantities"] = sp["quantity"].value_counts().to_dict()
    rep["sample_state"] = t["observations"]["sample_state"].value_counts().to_dict()
    rep["grids"] = len(dt.grids)
    rep["grid_ranges"] = sorted({(float(w[0]), float(w[-1]), len(w)) for w, _ in dt.grids.values()})[:5]
    rep["replicates_per_obs"] = sp.groupby(["observation_id", "quantity"]).size().describe()[
        ["mean", "max"]].round(2).to_dict()
    vals = np.concatenate([v[::25] for v in sp["values"]]) if len(sp) else np.array([np.nan])
    rep["value_quantiles"] = np.nanquantile(vals, [0.001, 0.5, 0.999]).round(4).tolist()
    tr = t["traits"]
    rep["traits"] = (tr.groupby("trait")["value"].describe()[["count", "min", "50%", "max"]]
                     .round(4).reset_index().to_dict("records") if len(tr) else [])
    rep["species_raw"] = int(dt.raw_species[[c for c in ("name", "code") if c in dt.raw_species]]
                             .drop_duplicates().shape[0])
    rep["treatments"] = t["treatments"][["factor", "level", "is_control"]].to_dict("records")[:20]
    rep["n_sites"], rep["n_plots"] = len(t["sites"]), len(t["plots"])
    return rep


def validate_partial(dt: DatasetTables) -> None:
    for name in ("datasets", "sites", "plots", "treatments", "entities", "observations",
                 "observation_treatments", "spectra", "traits"):
        schema.TABLES[name].validate(dt.tables[name], lazy=True)
    obs = set(dt.tables["observations"]["observation_id"])
    assert set(dt.tables["spectra"]["observation_id"]) <= obs
    tr = dt.tables["traits"]
    assert set(tr["observation_id"].dropna()) <= obs
    assert set(dt.tables["observations"]["entity_id"]) <= set(dt.tables["entities"]["entity_id"])


def write(dataset_id: str, dt: DatasetTables) -> None:
    d = interim_dir(dataset_id)
    for name, df in dt.tables.items():
        df.to_parquet(d / f"{name}.parquet", index=False)
    grids = pd.DataFrame({"grid_id": list(dt.grids),
                          "wavelength": [w for w, _ in dt.grids.values()],
                          "segment_start": [s for _, s in dt.grids.values()]})
    grids.to_parquet(d / "wavelength_grids.parquet", index=False)
    dt.raw_species.to_parquet(d / "raw_species.parquet", index=False)
    (d / "report.json").write_text(json.dumps(summarize(dt), indent=1, default=str))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["check", "dataset"])
    ap.add_argument("dataset_ids", nargs="+")
    a = ap.parse_args(argv)
    failed = []
    for ds in a.dataset_ids:
        try:
            dt = convert(ds)
            validate_partial(dt)
            if a.action == "dataset":
                write(ds, dt)
            print(json.dumps({ds: summarize(dt)}, indent=1, default=str))
        except Exception as e:  # report and continue with the others
            failed.append(ds)
            print(f"ERROR {ds}: {type(e).__name__}: {e}", file=sys.stderr)
    if failed:
        sys.exit(f"failed: {failed}")


if __name__ == "__main__":
    main()
