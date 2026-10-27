"""Combine per-batch fit outputs of all runs of a fit set.

    python -m leafrtm.aggregate --fitset-dir results/pilot

Writes ``<fitset>/summary/{summary,draws,spectra,info}.parquet`` with a
``run_id`` column and the run's model/errmodel/rt/data_use attributes.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from leafrtm.run import load_config


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fitset-dir", required=True)
    a = ap.parse_args(argv)
    d = Path(a.fitset_dir)
    cfg = load_config()
    runs = pd.DataFrame(cfg["fitsets"][d.name]["runs"]).rename(columns={"id": "run_id"})
    out = d / "summary"
    out.mkdir(parents=True, exist_ok=True)
    for part in ("summary", "draws", "spectra", "info"):
        frames = [pd.read_parquet(p) for p in sorted(d.glob(f"fits/*/batch*.{part}.parquet"))]
        frames = [f for f in frames if len(f)]
        df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        if len(df):
            df = df.merge(runs, on="run_id", how="left")
        df.to_parquet(out / f"{part}.parquet", index=False)
        print(part, len(df))


if __name__ == "__main__":
    main()
