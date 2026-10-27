"""Job files for the YCRC Dead Simple Queue (dSQ) on Bouchet (or any SLURM cluster).

    python -m leafrtm.hpc --fitset full --stage 1 [--pilot-dir results/pilot]

Writes ``hpc/fit_jobs_<fitset>_stage<k>.txt``: one line per (run, batch) whose
outputs are missing. Each line runs ``leafrtm.run`` for one batch; per-leaf
caching makes reruns cheap. With ``--pilot-dir`` it also prints the projected
core-hours from the pilot's median seconds per fit.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from leafrtm.run import ROOT, load_config


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fitset", default="full")
    ap.add_argument("--stage", type=int, choices=[1, 2], required=True)
    ap.add_argument("--python", default=".venv/bin/python")
    ap.add_argument("--pilot-dir", default=None)
    a = ap.parse_args(argv)
    cfg = load_config()
    fs = cfg["fitsets"][a.fitset]
    out = ROOT / "results" / a.fitset
    runs = [r for r in fs["runs"] if r["stage"] == a.stage]
    lines = []
    for r in runs:
        for b in range(fs["n_batches"]):
            if not (out / "fits" / r["id"] / f"batch{b:04d}.summary.parquet").exists():
                lines.append(f"cd {ROOT} && OMP_NUM_THREADS=1 {a.python} -W ignore -m leafrtm.run "
                             f"--fitset-dir results/{a.fitset} --run {r['id']} --batch {b}")
    hpc = ROOT / "hpc"
    hpc.mkdir(exist_ok=True)
    path = hpc / f"fit_jobs_{a.fitset}_stage{a.stage}.txt"
    path.write_text("\n".join(lines) + "\n")
    print(f"{len(lines)} jobs -> {path}")
    if a.pilot_dir and (out / "fit_set.parquet").exists():
        info = pd.read_parquet(Path(a.pilot_dir) / "summary" / "info.parquet")
        sec = info[info["status"] == "ok"].groupby("model")["seconds"].median()
        fit_set = pd.read_parquet(out / "fit_set.parquet")
        n_rt = int(fit_set["has_RT"].sum())
        hours = 0.0
        for r in runs:
            n = n_rt if r.get("only_rt") else len(fit_set)
            hours += n * sec.get(r["model"], sec.median()) / 3600
        print(f"projected ~{hours:,.0f} core-hours for stage {a.stage} "
              f"(pilot medians per model; sampler settings may differ)")


if __name__ == "__main__":
    main()
