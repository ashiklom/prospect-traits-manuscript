"""Learn error-model structure from first-stage (``ou``) fits.

    python -m leafrtm.learn --fitset-dir results/pilot --stage1-runs <id> [<id> ...]

Writes ``<fitset>/discrepancy/ss_correlation.json`` (Matern-3/2 length scale
and nugget per model and quantity, from residuals standardized by each leaf's
posterior heteroskedastic scale) and ``<fitset>/discrepancy/<model>/<tag>/
<class>/fold<k>.npz`` (learned discrepancies; tag ``rt`` = joint R/T factors,
``sep`` = separate R and T factors; class = sample state x foreoptic; fold =
dataset-level cross-fitting fold: leaves of fold k are learned from).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from leafrtm import discrepancy as disc
from leafrtm.prep import read_batch
from leafrtm.run import discrepancy_class, fold_of, load_config

MIN_LEAVES = 15


def residuals(fitset_dir: Path, run_id: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    d = Path(fitset_dir) / "fits" / run_id
    spec = pd.concat([pd.read_parquet(p) for p in sorted(d.glob("*.spectra.parquet"))],
                     ignore_index=True)
    summ = pd.concat([pd.read_parquet(p) for p in sorted(d.glob("*.summary.parquet"))],
                     ignore_index=True)
    spec = spec[spec["fitted"]].copy()
    spec["resid"] = spec["y"] - spec["yhat_mean"]
    return spec, summ


def leaf_meta(fitset_dir: Path) -> pd.DataFrame:
    rows = []
    for p in sorted((Path(fitset_dir) / "inputs").glob("*.pkl")):
        b = read_batch(p)
        m = b["meta"][["observation_id", "dataset_id", "sample_state"]].copy()
        m["foreoptic"] = [b["blocks"][o][0].instrument_class if b["blocks"].get(o) else "unknown"
                          for o in m["observation_id"]]
        rows.append(m)
    m = pd.concat(rows, ignore_index=True)
    m["klass"] = [discrepancy_class(s, f) for s, f in zip(m["sample_state"], m["foreoptic"])]
    return m


def standardized(spec, summ, model):
    """Residual / (sigma0 + sigma1 |yhat|) using posterior means per quantity."""
    sig = summ[summ["param"].str.match(r"sigma[01]\[\d\]")].pivot_table(
        index="observation_id", columns="param", values="mean")
    out = []
    for (oid, q), g in spec.groupby(["observation_id", "quantity"]):
        qs = sorted(spec.loc[spec.observation_id == oid, "quantity"].unique(),
                    key=["R", "T", "A", "PA"].index)
        j = qs.index(q)
        s0 = sig.at[oid, f"sigma0[{j}]"] if f"sigma0[{j}]" in sig else np.nan
        s1 = sig.at[oid, f"sigma1[{j}]"] if f"sigma1[{j}]" in sig else 0.0
        scale = s0 + s1 * np.abs(g["yhat_mean"].to_numpy())
        out.append((q, g["wavelength"].to_numpy(), g["resid"].to_numpy() / scale))
    return out


def regular_series(wl, z, step=1.0):
    """Residual series on its own regular sub-grid (for the correlogram)."""
    grid = np.arange(np.ceil(wl.min()), np.floor(wl.max()) + step, step)
    v = np.interp(grid, wl, z)
    gap = np.diff(wl).max() if len(wl) > 1 else 0
    return (grid, v) if gap <= 5 * step else (wl, z)


def learn_all(fitset_dir: Path, stage1_runs: list[str], k: int, n_folds: int) -> None:
    cfg_out = Path(fitset_dir) / "discrepancy"
    cfg_out.mkdir(parents=True, exist_ok=True)
    meta = leaf_meta(fitset_dir).set_index("observation_id")
    corr = {}
    report = {}
    for run_id in stage1_runs:
        spec, summ = residuals(fitset_dir, run_id)
        model = summ["run_id"].iloc[0].split("-")[0] if "-" in run_id else run_id
        model = run_id.split("_")[0]
        fresh = spec["observation_id"].map(meta["sample_state"]) == "fresh"
        series = standardized(spec[fresh], summ, model)
        for q in ("R", "T", "all"):
            ss = [regular_series(w, z) for qq, w, z in series if q in ("all", qq)]
            if len(ss) >= 5:
                ell, nug, _, _ = disc.fit_correlogram(ss)
                corr[f"{model}/{q}"] = {"ell": ell, "nugget": nug, "n_series": len(ss)}
        leaves = {}
        for oid, g in spec.groupby("observation_id"):
            leaves[oid] = {q: (h["wavelength"].to_numpy(), h["resid"].to_numpy())
                           for q, h in g.groupby("quantity") if q in ("R", "T")}
        for tag in ("rt", "sep"):
            for klass in sorted(set(meta["klass"])) + ["all"]:
                for fold in range(n_folds):
                    ids = [o for o in leaves if fold_of(o, n_folds) != fold
                           and (klass == "all" or meta.at[o, "klass"] == klass)]
                    if len(ids) < MIN_LEAVES:
                        continue
                    try:
                        d = (disc.learn([leaves[o] for o in ids], k=k) if tag == "rt"
                             else learn_separate([leaves[o] for o in ids], k))
                    except ValueError:
                        continue
                    p = cfg_out / model / tag / klass / f"fold{fold}.npz"
                    p.parent.mkdir(parents=True, exist_ok=True)
                    d.save(p)
                    report[f"{model}/{tag}/{klass}/fold{fold}"] = {
                        "n_leaves": len(ids), "ell": d.ell, "nugget": d.nugget}
    (cfg_out / "ss_correlation.json").write_text(json.dumps(corr, indent=1))
    (cfg_out / "report.json").write_text(json.dumps(report, indent=1))


def learn_separate(leaves, k) -> disc.Discrepancy:
    """R and T learned independently, combined block-diagonally."""
    parts = []
    for q in ("R", "T"):
        sub = [{q: lf[q]} for lf in leaves if q in lf]
        parts.append(disc.learn(sub, k=k, quantities=(q,)) if len(sub) >= k + 2 else None)
    ng = len(disc.REF_GRID)
    kk = sum(p.Lambda.shape[1] for p in parts if p is not None)
    mu, D, Lam = np.zeros(2 * ng), np.ones(2 * ng), np.zeros((2 * ng, kk))
    col = 0
    ells, nugs, n = [], [], 0
    for i, p in enumerate(parts):
        if p is None:
            D[i * ng:(i + 1) * ng] = np.nan
            continue
        sl = slice(i * ng, (i + 1) * ng)
        mu[sl], D[sl] = p.mu, p.D
        Lam[sl, col:col + p.Lambda.shape[1]] = p.Lambda
        col += p.Lambda.shape[1]
        ells.append(p.ell)
        nugs.append(p.nugget)
        n = max(n, p.n_leaves)
    if np.isnan(D).any():  # quantity without enough leaves: borrow the other's scale
        D = np.where(np.isnan(D), np.nanmedian(D), D)
    return disc.Discrepancy(("R", "T"), disc.REF_GRID, mu, Lam, D, float(np.mean(ells)),
                            float(np.mean(nugs)), n)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fitset-dir", required=True)
    ap.add_argument("--stage1-runs", nargs="+", required=True)
    a = ap.parse_args(argv)
    cfg = load_config()
    learn_all(Path(a.fitset_dir), a.stage1_runs, cfg["discrepancy"]["k"],
              cfg["discrepancy"]["n_folds"])


if __name__ == "__main__":
    main()
