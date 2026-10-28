"""Select the observations of a fit set and write per-batch inversion inputs.

    python -m leafrtm.select --fitset pilot
"""

from __future__ import annotations

import argparse
from pathlib import Path

from leafrtm.db.harmonize import load
from leafrtm.prep import blocks_for, observation_table, select, write_batch
from leafrtm.run import ROOT, load_config


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fitset", required=True)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    cfg = load_config()
    fs = cfg["fitsets"][a.fitset]
    out = Path(a.out or ROOT / "results" / a.fitset)
    db = load()
    obs = observation_table(db, quantities=fs.get("quantities", ["R", "T"]))
    chosen = select(obs, fs)
    out.mkdir(parents=True, exist_ok=True)
    chosen.to_parquet(out / "fit_set.parquet", index=False)
    blocks = blocks_for(chosen["observation_id"], db, mask_ranges=cfg.get("mask_ranges", []),
                        quantities=fs.get("quantities", ["R", "T"]))
    for b in range(fs["n_batches"]):
        m = chosen[chosen["batch"] == b]
        write_batch(out / "inputs" / f"{b:04d}.pkl", {o: blocks[o] for o in m["observation_id"]
                                                  if o in blocks}, m)
    print(f"{len(chosen)} observations in {fs['n_batches']} batches -> {out}")


if __name__ == "__main__":
    main()
