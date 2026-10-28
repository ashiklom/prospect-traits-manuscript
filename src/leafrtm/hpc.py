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
    if a.pilot_dir:
        print(f"projected ~{projected_hours(cfg, a.fitset, runs, Path(a.pilot_dir)):,.0f} core-hours "
              f"for stage {a.stage}")


def projected_hours(cfg: dict, fitset: str, runs: list[dict], pilot_dir: Path) -> float:
    """Core-hours from the pilot's mean seconds per fit and model. The sampling part
    is scaled by the ratio of gradient-evaluation budgets (chains x (tune + draws)) of
    the two fit sets; the mode search and compilation are not. PROSPECT versions not
    in the pilot use PROSPECT-D's timings."""
    from leafrtm.db.harmonize import load
    from leafrtm.prep import observation_table, select

    fs, pilot = cfg["fitsets"][fitset], cfg["fitsets"][pilot_dir.name]
    budget = lambda s: s["chains"] * (s["tune"] + s["draws"])  # noqa: E731
    ratio = budget(fs["sampler"]) / budget(pilot["sampler"])
    info = pd.read_parquet(pilot_dir / "summary" / "info.parquet")
    info = info[info["status"] == "ok"]
    per_fit = ((info["seconds"] - info["sampling_seconds"]) + ratio * info["sampling_seconds"])
    sec = per_fit.groupby(info["model"]).mean()
    chosen = select(observation_table(load(), quantities=fs.get("quantities", ["R", "T"])), fs)
    n_all, n_rt = len(chosen), int(chosen["has_RT"].sum())
    hours = sum((n_rt if r.get("only_rt") else n_all)
                * sec.get(r["model"], sec.get("prospectd", sec.mean())) / 3600 for r in runs)
    print(f"{n_all:,} observations ({n_rt:,} with R and T); sampler budget x{ratio:.2f} of the pilot's")
    return hours


if __name__ == "__main__":
    main()
