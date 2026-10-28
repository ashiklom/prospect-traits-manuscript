"""Fit one batch of observations for one run specification.

    python -m leafrtm.run --fitset-dir results/pilot --run <run_id> --batch 3

A *run* (``config/config.yaml``) fixes the RTM, error model, R/T correlation,
data use ("all" quantities or "R" only) and stage. Per-observation results
are cached under ``<fitset>/cache/<run_id>/`` (resumable); the batch output
is three parquet files: posterior summaries, thinned draws, fitted spectra.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import time
import traceback
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import yaml

from leafrtm import discrepancy as disc
from leafrtm import likelihood as L
from leafrtm.fit import FitSpec, Inverter, data_arrays, load_priors
from leafrtm.obsop import QCODE, make_obsdata
from leafrtm.prep import data_use_variants, read_batch
from leafrtm.rtm import get_model

ROOT = Path(__file__).resolve().parents[2]
N_KEEP_DRAWS = 100


def load_config(path=ROOT / "config" / "config.yaml") -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def discrepancy_class(sample_state: str, foreoptic: str) -> str:
    """Leaves share a learned discrepancy within (fresh or not) x foreoptic group."""
    optic = {"integrating_sphere": "sphere", "leaf_clip": "clip", "contact_probe": "clip"}.get(
        foreoptic, "other")
    return f"{'fresh' if sample_state == 'fresh' else 'dry'}_{optic}"


def fold_of(observation_id: str, n_folds: int) -> int:
    """Cross-fitting fold by dataset (all leaves of a dataset share a fold)."""
    ds = observation_id.split(":")[0]
    return int(hashlib.md5(ds.encode()).hexdigest(), 16) % n_folds


# ---------------------------------------------------------------------------
# Error-model inputs
# ---------------------------------------------------------------------------


def ss_arrays(od, corr: dict, model: str):
    """Vecchia factors per block with the learned Matern (ell, nugget) for
    (model, quantity); falls back to the model's pooled values."""
    m = 30
    nb = np.zeros((len(od), m), dtype=np.int64)
    b = np.zeros((len(od), m))
    S = np.ones(len(od))
    for bi in np.unique(od.block):
        idx = np.flatnonzero(od.block == bi)
        q = [k for k, v in QCODE.items() if v == od.qcode[idx[0]]][0]
        p = corr.get(f"{model}/{q}") or corr[f"{model}/all"]
        n_, b_, S_ = L.vecchia(od.wavelength[idx], od.start[idx],
                               lambda d, p=p: L.matern32(d, p["ell"], p["nugget"]), m=m)
        nb[idx], b[idx], S[idx] = idx[n_], b_, S_
    return {"nb": nb, "b": b, "S": S}


def ld_arrays(od, klass, discrepancies: dict, model: str, oid: str, n_folds: int, rt: bool):
    fold = fold_of(oid, n_folds)
    tag = "rt" if rt else "sep"
    for key in (f"{model}/{tag}/{klass}/fold{fold}", f"{model}/{tag}/all/fold{fold}"):
        if key in discrepancies:
            return disc.map_to_obs(discrepancies[key], od), key
    raise KeyError(f"No discrepancy for {model}/{tag}/{klass}/fold{fold}")


# ---------------------------------------------------------------------------
# Posterior summaries and predictions
# ---------------------------------------------------------------------------


def _diag(trace, fn, **kw):
    d = fn(trace, **kw)
    return d.posterior if hasattr(d, "posterior") else d


def summarize(trace, oid, run_id, rtm_params) -> pd.DataFrame:
    post = trace.posterior
    rows = []
    ess_bulk = _diag(trace, az.ess, method="bulk")
    ess_tail = _diag(trace, az.ess, method="tail")
    rhat = _diag(trace, az.rhat)

    def val(ds, name, idx):
        x = np.asarray(ds[name])
        return float(x[idx] if x.ndim else x)

    for name in post.data_vars:
        x = np.asarray(post[name])  # chain, draw, ...
        flat = x.reshape(x.shape[0] * x.shape[1], -1)
        for j in range(flat.shape[1]):
            idx = np.unravel_index(j, x.shape[2:]) if x.ndim > 2 else ()
            v = flat[:, j]
            q = np.quantile(v, [0.05, 0.25, 0.5, 0.75, 0.95])
            rows.append({
                "observation_id": oid, "run_id": run_id,
                "param": name + (f"[{j}]" if x.ndim > 2 else ""),
                "is_rtm": name in rtm_params,
                "mean": v.mean(), "sd": v.std(), "q05": q[0], "q25": q[1], "q50": q[2],
                "q75": q[3], "q95": q[4],
                "ess_bulk": val(ess_bulk, name, idx), "ess_tail": val(ess_tail, name, idx),
                "rhat": val(rhat, name, idx)})
    return pd.DataFrame(rows)


def boundary_stats(model, draws: pd.DataFrame, priors: dict) -> dict:
    """Where the posterior sits relative to hard limits.

    ``wall_margin`` / ``wall_margin_q05``: the model's validity margin (e.g. the
    smallest LIBERTY discriminant) at the posterior median and its 5 % quantile
    over the thinned draws. ``frac_at_bound``: share of draws with any RTM
    parameter within 1 % (of its range) of a hard prior bound."""
    out = {}
    med = {p: float(draws[p].median()) for p in model.params}
    if model.validity_margin(med) is not None:
        m = np.array([model.validity_margin({p: r[p] for p in model.params})
                      for _, r in draws.iterrows()])
        out["wall_margin"] = model.validity_margin(med)
        out["wall_margin_q05"] = float(np.quantile(m, 0.05))
    near = np.zeros(len(draws), bool)
    for p in model.params:
        spec = priors["rtm"][p]
        lo, hi = spec.get("lower"), spec.get("upper")
        if lo is None and hi is None:
            continue
        x = draws[p].to_numpy()
        tol = 0.01 * ((hi - lo) if lo is not None and hi is not None else abs(lo if hi is None else hi))
        if lo is not None:
            near |= x - lo < tol
        if hi is not None:
            near |= hi - x < tol
    out["frac_at_bound"] = float(near.mean())
    return out


def thinned_draws(trace, params, n=N_KEEP_DRAWS, seed=0) -> pd.DataFrame:
    post = trace.posterior
    total = post.sizes["chain"] * post.sizes["draw"]
    idx = np.random.default_rng(seed).choice(total, size=min(n, total), replace=False)
    out = {"draw": np.arange(len(idx))}
    for p in params:
        out[p] = np.asarray(post[p]).reshape(-1)[idx]
    return pd.DataFrame(out)


def predict(model, draws: pd.DataFrame, od) -> np.ndarray:
    """(n_draws, n_bands) predicted spectra on the observation's bands."""
    out = np.empty((len(draws), len(od)))
    for i, row in enumerate(draws[list(model.params)].itertuples(index=False)):
        R, T = model.simulate(dict(zip(model.params, row)))
        out[i] = od.predict_np(R, T)
    return out


def spectra_summary(model, draws, od_all, fitted_q, oid, run_id) -> pd.DataFrame:
    pred = predict(model, draws, od_all)
    qnames = {v: k for k, v in QCODE.items()}
    q = np.nanquantile(pred, [0.05, 0.5, 0.95], axis=0)
    return pd.DataFrame({
        "observation_id": oid, "run_id": run_id, "block": od_all.block,
        "quantity": [qnames[c] for c in od_all.qcode],
        "wavelength": od_all.wavelength, "y": od_all.y,
        "yhat_mean": np.nanmean(pred, axis=0), "yhat_q05": q[0], "yhat_q50": q[1],
        "yhat_q95": q[2],
        "fitted": [qnames[c] in fitted_q for c in od_all.qcode]})


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------


def run_batch(fitset_dir: Path, run: dict, batch: int, cfg: dict) -> dict:
    fitset_dir = Path(fitset_dir)
    run_id = run["id"]
    sampler = {**cfg["sampler"], **run.get("sampler", {})}
    inputs = read_batch(fitset_dir / "inputs" / f"{batch:04d}.pkl")
    model = get_model(run["model"])
    cache = fitset_dir / "cache" / run_id
    cache.mkdir(parents=True, exist_ok=True)
    corr, discs = {}, {}
    if run["errmodel"] == "ss":
        corr = json.loads((fitset_dir / "discrepancy" / "ss_correlation.json").read_text())
    if run["errmodel"] == "ld":
        for p in (fitset_dir / "discrepancy").glob(f"{run['model']}/*/*/fold*.npz"):
            rel = p.relative_to(fitset_dir / "discrepancy").with_suffix("")
            discs[str(rel)] = disc.Discrepancy.load(p)
    priors = load_priors()
    inverters: dict[str, Inverter] = {}
    status = []
    meta = inputs["meta"].set_index("observation_id")
    for oid in meta.index:
        fn = cache / (hashlib.md5(oid.encode()).hexdigest() + ".pkl")
        if fn.exists():
            continue
        all_blocks = inputs["blocks"].get(oid, [])
        if run.get("only_rt") and not meta.at[oid, "has_RT"]:
            continue
        variants = data_use_variants(all_blocks, [run["data_use"]])
        if run["data_use"] not in variants:
            continue
        blocks = [b for b in variants[run["data_use"]]
                  if b.quantity in cfg.get("fit_quantities", ["R", "T", "A", "PA"])]
        if not blocks:
            continue
        t0 = time.perf_counter()
        try:
            od = make_obsdata(blocks, model.wavelength)
            od_all = make_obsdata(all_blocks, model.wavelength)
            spec = FitSpec(run["model"], run["errmodel"], run.get("rt", False), od.quantities)
            extra, dkey = None, None
            if spec.errmodel == "ss":
                extra = ss_arrays(od, corr, run["model"])
            elif spec.errmodel == "ld":
                klass = discrepancy_class(meta.at[oid, "sample_state"], blocks[0].instrument_class)
                extra, dkey = ld_arrays(od, klass, discs, run["model"], oid,
                                        cfg["discrepancy"]["n_folds"], spec.rt)
            arrays = data_arrays(od, spec, extra)
            if spec.key not in inverters:
                inverters[spec.key] = Inverter(spec, arrays)
            seed = zlib_seed(oid, run_id)
            trace, t_sample = inverters[spec.key].sample(
                arrays, draws=sampler["draws"], tune=sampler["tune"], chains=sampler["chains"],
                seed=seed, target_accept=sampler.get("target_accept", 0.8))
            summ = summarize(trace, oid, run_id, model.params)
            div = int(np.asarray(trace.sample_stats["diverging"]).sum())
            draws = thinned_draws(trace, list(model.params), seed=seed)
            draws.insert(0, "observation_id", oid)
            draws.insert(1, "run_id", run_id)
            spec_sum = spectra_summary(model, draws, od_all, set(od.quantities), oid, run_id)
            n_grad = int(np.asarray(trace.sample_stats["n_steps"]).sum()) \
                if "n_steps" in trace.sample_stats else None
            info = {"observation_id": oid, "run_id": run_id, "status": "ok",
                    "seconds": time.perf_counter() - t0, "sampling_seconds": t_sample,
                    "n_grad": n_grad, "divergences": div, "n_bands": len(od),
                    "quantities": "".join(od.quantities), "discrepancy": dkey,
                    "max_rhat": float(summ.loc[summ.is_rtm, "rhat"].max()),
                    "min_ess_bulk": float(summ.loc[summ.is_rtm, "ess_bulk"].min()),
                    "min_ess_tail": float(summ.loc[summ.is_rtm, "ess_tail"].min()),
                    **boundary_stats(model, draws, priors)}
            payload = {"summary": summ, "draws": draws, "spectra": spec_sum, "info": info}
        except Exception as e:  # record and continue
            info = {"observation_id": oid, "run_id": run_id, "status": f"error: {e}",
                    "seconds": time.perf_counter() - t0, "traceback": traceback.format_exc()}
            payload = {"info": info}
        with open(fn, "wb") as f:
            pickle.dump(payload, f)
        status.append(info)
    return {"n": len(status)}


def zlib_seed(oid, run_id) -> int:
    return int(hashlib.md5(f"{oid}|{run_id}".encode()).hexdigest()[:8], 16)


def collect(fitset_dir: Path, run_id: str, batch_obs: list[str], out_prefix: Path) -> None:
    """Gather cached per-observation results of one batch into parquet files."""
    cache = Path(fitset_dir) / "cache" / run_id
    parts = {"summary": [], "draws": [], "spectra": [], "info": []}
    for oid in batch_obs:
        fn = cache / (hashlib.md5(oid.encode()).hexdigest() + ".pkl")
        if not fn.exists():
            continue
        with open(fn, "rb") as f:
            p = pickle.load(f)
        for k in parts:
            if k in p:
                parts[k].append(p[k] if k != "info" else pd.DataFrame([p[k]]))
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    for k, v in parts.items():
        df = pd.concat(v, ignore_index=True) if v else pd.DataFrame()
        df.to_parquet(f"{out_prefix}.{k}.parquet", index=False)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fitset-dir", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--batch", type=int, required=True)
    ap.add_argument("--config", default=str(ROOT / "config" / "config.yaml"))
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    fitset = Path(a.fitset_dir).name
    run = next(r for r in cfg["fitsets"][fitset]["runs"] if r["id"] == a.run)
    cfg_fs = {**cfg, **cfg["fitsets"][fitset]}
    run_batch(Path(a.fitset_dir), run, a.batch, cfg_fs)
    meta = read_batch(Path(a.fitset_dir) / "inputs" / f"{a.batch:04d}.pkl")["meta"]
    collect(Path(a.fitset_dir), a.run, list(meta["observation_id"]),
            Path(a.fitset_dir) / "fits" / a.run / f"batch{a.batch:04d}")


if __name__ == "__main__":
    main()
