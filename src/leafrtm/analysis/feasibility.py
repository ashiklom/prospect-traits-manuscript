"""Can absorptance (A = 1 - R - T) and pseudo-absorbance (PA = log10 1/R) be
fitted? Compares fits of A alone with fits of R + T on the same leaves
(pepper), and reports convergence and plausibility of PA fits (ACCP).

    python -m leafrtm.analysis.feasibility --fitset-dir results/feasibility --out results/pilot/analysis
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from leafrtm.analysis.estimates import measured_traits, summarize_draws, trait_draws
from leafrtm.db.harmonize import load


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fitset-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    d = Path(a.fitset_dir)
    draws = pd.read_parquet(d / "summary" / "draws.parquet")
    info = pd.read_parquet(d / "summary" / "info.parquet")
    est = summarize_draws(trait_draws(draws))
    meas = measured_traits(load())
    est = est.merge(meas, on=["observation_id", "trait"], how="left")
    ok = info[info["status"] == "ok"]
    conv = ok.assign(converged=(ok["max_rhat"] < 1.01) & (ok["min_ess_bulk"] >= 200)).groupby(
        "run_id").agg(n=("observation_id", "size"), converged=("converged", "mean"),
                      median_seconds=("seconds", "median"))
    rows = []
    a_ = est[est["run_id"] == "prospectd_ou_A"]
    rt = est[est["run_id"] == "prospectd_ou_RT"]
    j = a_.merge(rt, on=["observation_id", "trait"], suffixes=("_A", "_RT"))
    for t, g in j.groupby("trait"):
        ok_ = (g["q50_A"] > 0) & (g["q50_RT"] > 0)
        if ok_.sum() < 3:
            continue
        rows.append({"comparison": "A vs R+T (same leaves)", "trait": t, "n": int(ok_.sum()),
                     "r2_log": float(np.corrcoef(np.log(g.loc[ok_, "q50_A"]), np.log(g.loc[ok_, "q50_RT"]))[0, 1] ** 2),
                     "median_ratio": float(np.median(g.loc[ok_, "q50_A"] / g.loc[ok_, "q50_RT"])),
                     "intervals_overlap": float(np.mean((g["q05_A"] <= g["q95_RT"]) & (g["q05_RT"] <= g["q95_A"])))})
    pa = est[est["run_id"] == "prospectd_ou_PA"].dropna(subset=["measured"])
    for t, g in pa.groupby("trait"):
        ok_ = (g["q50"] > 0) & (g["measured"] > 0)
        if ok_.sum() < 3:
            continue
        rows.append({"comparison": "PA vs measured", "trait": t, "n": int(ok_.sum()),
                     "r2_log": float(np.corrcoef(np.log(g.loc[ok_, "q50"]), np.log(g.loc[ok_, "measured"]))[0, 1] ** 2),
                     "median_ratio": float(np.median(g.loc[ok_, "q50"] / g.loc[ok_, "measured"])),
                     "intervals_overlap": float(np.mean((g["q05"] <= g["measured"]) & (g["measured"] <= g["q95"])))})
    out = pd.DataFrame(rows)
    for run_id, r in conv.iterrows():
        out = pd.concat([out, pd.DataFrame([{"comparison": f"convergence {run_id}", "n": r["n"],
                                             "converged_share": r["converged"],
                                             "median_seconds": r["median_seconds"]}])], ignore_index=True)
    Path(a.out).mkdir(parents=True, exist_ok=True)
    out.to_parquet(Path(a.out) / "apa_feasibility.parquet", index=False)
    print(out)


if __name__ == "__main__":
    main()
