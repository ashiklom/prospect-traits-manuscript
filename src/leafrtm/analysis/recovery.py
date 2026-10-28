"""Synthetic recovery: simulate spectra from each RTM with known parameters and
correlated heteroskedastic noise, fit them, and check bias and calibration.

    python -m leafrtm.analysis.recovery --out results/pilot/analysis [--n 8]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from leafrtm import likelihood as L
from leafrtm.fit import FitSpec, Inverter, data_arrays, load_priors
from leafrtm.obsop import Block, make_obsdata
from leafrtm.rtm import get_model

TRUE_ELL, TRUE_NUGGET = 40.0, 0.05
SIGMA0, SIGMA1 = 0.003, 0.02


def draw_theta(model, priors, rng):
    """Parameters from the central 10-90 % of each prior."""
    from scipy import stats
    out = {}
    for p in model.params:
        s = priors["rtm"][p]
        u = rng.uniform(0.1, 0.9)
        if s["dist"] == "lognormal":
            out[p] = float(np.exp(s["meanlog"] + s["sdlog"] * stats.norm.ppf(u)))
        elif s["dist"] == "truncnormal":
            lo = (s.get("lower", -np.inf) - s["mu"]) / s["sigma"]
            hi = (s.get("upper", np.inf) - s["mu"]) / s["sigma"]
            out[p] = float(stats.truncnorm(lo, hi, s["mu"], s["sigma"]).ppf(u))
        elif s["dist"] == "uniform":
            out[p] = float(s["lower"] + u * (s["upper"] - s["lower"]))
    return out


def simulate(model, theta, rng):
    R, T = model.simulate(theta)
    wl = model.wavelength
    keep = np.isfinite(R) & np.isfinite(T)
    blocks = []
    corr = lambda d: L.matern32(d, TRUE_ELL, TRUE_NUGGET)  # noqa: E731
    C = L.dense_corr(wl[keep], np.r_[True, np.zeros(keep.sum() - 1, bool)], corr)
    Lc = np.linalg.cholesky(C + 1e-9 * np.eye(len(C)))
    for q, y in (("R", R), ("T", T)):
        s = SIGMA0 + SIGMA1 * y[keep]
        blocks.append(Block(q, wl[keep], y[keep] + s * (Lc @ rng.normal(size=len(C)))))
    return blocks


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--models", nargs="+", default=["prospectd", "liberty", "leafsip"])
    a = ap.parse_args(argv)
    priors = load_priors()
    rng = np.random.default_rng(42)
    rows = []
    for name in a.models:
        model = get_model(name)
        inv = None
        for i in range(a.n):
            theta = draw_theta(model, priors, rng)
            blocks = simulate(model, theta, rng)
            od = make_obsdata(blocks, model.wavelength)
            spec = FitSpec(name, od.quantities)
            arr = data_arrays(od, spec)
            inv = inv or Inverter(spec, arr)
            tr, secs = inv.sample(arr, draws=500, tune=500, chains=4, seed=i)
            for p in model.params:
                x = np.asarray(tr.posterior[p]).ravel()
                q05, q50, q95 = np.quantile(x, [0.05, 0.5, 0.95])
                rows.append({"model": name, "errmodel": "ou", "replicate": i, "param": p,
                             "true": theta[p], "q05": q05, "q50": q50, "q95": q95,
                             "covered90": q05 <= theta[p] <= q95,
                             "log_error": float(np.log(q50 / theta[p])) if theta[p] > 0 else np.nan,
                             "seconds": secs})
    df = pd.DataFrame(rows)
    summ = df.groupby(["model", "errmodel", "param"]).agg(
        coverage90=("covered90", "mean"), median_log_error=("log_error", "median"),
        rmse_log=("log_error", lambda x: float(np.sqrt(np.nanmean(x**2)))), n=("replicate", "nunique"),
        median_seconds=("seconds", "median")).reset_index()
    Path(a.out).mkdir(parents=True, exist_ok=True)
    df.to_parquet(Path(a.out) / "recovery_draws.parquet", index=False)
    summ.to_parquet(Path(a.out) / "recovery.parquet", index=False)
    print(summ)


if __name__ == "__main__":
    main()
