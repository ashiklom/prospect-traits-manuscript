"""All analyses, from aggregated fits + the harmonized database.

    python -m leafrtm.analysis.compute --fitset-dir results/pilot

Writes tidy parquet tables to ``<fitset>/analysis/`` that the website reads.
Part A (RTM evaluation): intercompare, validation_leaf, validation_species,
spectral_residuals, heldout, error_model_comparison, diagnostics,
identifiability, sample_state, protocol. Part B (trait variation; uses the
best-estimate product): best_config, best_estimates, variance, treatments,
correlations, within_species_correlations.
"""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from leafrtm.analysis.estimates import (crps_samples, measured_traits, plant_traits,
                                        summarize_draws, trait_draws)
from leafrtm.db.harmonize import load
from leafrtm.fit import load_priors

MIN_N = 5


def with_runs(df: pd.DataFrame, runs: pd.DataFrame) -> pd.DataFrame:
    """Attach run attributes (model, errmodel, rt, data_use) not already present."""
    cols = ["run_id"] + [c for c in runs.columns if c != "run_id" and c not in df.columns]
    return df.merge(runs[cols], on="run_id", how="left")


# ---------------------------------------------------------------------------
# Context
# ---------------------------------------------------------------------------


def observation_context(db: dict) -> pd.DataFrame:
    obs = db["observations"].merge(db["entities"][["entity_id", "species_id", "plant_id", "plot_id"]],
                                   on="entity_id", how="left")
    sp = db["spectra"].groupby("observation_id")["foreoptic"].agg(lambda s: s.mode().iat[0])
    obs = obs.merge(sp.rename("foreoptic"), on="observation_id", how="left")
    obs = obs.merge(db["species"][["species_id", "scientific_name", "family", "growth_form",
                                   "leaf_type"]], on="species_id", how="left")
    obs["year"] = pd.to_datetime(obs["date"], errors="coerce").dt.year.astype("Int64")
    return obs


# ---------------------------------------------------------------------------
# Suspect groups and exclusions
# ---------------------------------------------------------------------------

FLAG_DEFAULTS = {"min_n": 3, "max_abs_log_deviation": 0.405, "min_share_runs": 0.5}


def suspect_groups(leaf: pd.DataFrame, ctx: pd.DataFrame, flag: dict | None = None) -> pd.DataFrame:
    """(dataset, year, trait) groups whose RTM estimates miss the measurements far more
    than the same trait does elsewhere.

    Per run and trait, a group's median log(estimate / measured) is compared with the
    median over all *other* groups of that run and trait, so biases common to a trait
    everywhere (an RTM problem) are not blamed on one dataset. A group is flagged when it
    has at least ``min_n`` leaves and its deviation exceeds ``max_abs_log_deviation``
    (0.405 = a factor of 1.5) in at least ``min_share_runs`` of the runs estimating the
    trait. Consistent misses across many related spectra point to the measurements (units,
    method) or the spectra of that group rather than to individual leaves."""
    f = {**FLAG_DEFAULTS, **(flag or {})}
    v = leaf[~leaf["relative"] & (leaf["q50"] > 0) & (leaf["measured"] > 0)].merge(
        ctx[["observation_id", "year"]], on="observation_id", how="left")
    v = v.assign(lr=np.log(v["q50"] / v["measured"]), year=v["year"].fillna(-1))
    g = (v.groupby(["run_id", "trait", "dataset_id", "year"])
         .agg(n=("lr", "size"), med_log_ratio=("lr", "median"), coverage90=("covered90", "mean"))
         .reset_index())
    rows = []
    for (run, trait), h in g.groupby(["run_id", "trait"]):
        for i, r in h.iterrows():
            others = h.drop(index=i)
            base = others.loc[others["n"] >= f["min_n"], "med_log_ratio"].median() if len(others) else np.nan
            rows.append({**r.to_dict(), "baseline_log_ratio": base, "deviation": r["med_log_ratio"] - base})
    g = pd.DataFrame(rows)
    if g.empty:
        return g
    g["flag_run"] = (g["n"] >= f["min_n"]) & (g["deviation"].abs() > f["max_abs_log_deviation"])
    out = (g.groupby(["dataset_id", "year", "trait"])
           .agg(n_leaves=("n", "max"), n_runs=("run_id", "nunique"), share_runs_flagged=("flag_run", "mean"),
                median_ratio=("med_log_ratio", lambda x: float(np.exp(np.median(x)))),
                deviation_ratio=("deviation", lambda x: float(np.exp(np.nanmedian(x)))),
                coverage90=("coverage90", "median"))
           .reset_index())
    out["flagged"] = (out["n_leaves"] >= f["min_n"]) & (out["share_runs_flagged"] >= f["min_share_runs"])
    out["year"] = out["year"].astype(int).replace(-1, pd.NA)
    return out.sort_values(["flagged", "share_runs_flagged"], ascending=False).reset_index(drop=True)


def exclusion_rules(analysis_cfg: dict, suspects: pd.DataFrame) -> pd.DataFrame:
    """Exclusions from the config (``analysis.exclude``) plus, with
    ``analysis.auto_exclude_flagged``, every flagged suspect group (scope "ecology").
    Columns: dataset_id, year (NA = all years), trait (NA = all traits), scope, reason, source."""
    rows = []
    for e in analysis_cfg.get("exclude", []) or []:
        for t in e.get("traits") or [None]:
            rows.append({"dataset_id": e["dataset"], "year": e.get("year"), "trait": t,
                         "scope": e.get("scope", "ecology"), "reason": e.get("reason", ""), "source": "config"})
    if analysis_cfg.get("auto_exclude_flagged") and len(suspects):
        for r in suspects[suspects["flagged"]].itertuples():
            rows.append({"dataset_id": r.dataset_id, "year": r.year, "trait": r.trait, "scope": "ecology",
                         "reason": f"auto: estimates/measured {r.deviation_ratio:.2f}x the trait's ratio elsewhere",
                         "source": "flag"})
    cols = ["dataset_id", "year", "trait", "scope", "reason", "source"]
    return pd.DataFrame(rows, columns=cols)


def excluded(df: pd.DataFrame, ctx: pd.DataFrame, rules: pd.DataFrame, scopes=("ecology", "all")) -> pd.Series:
    """Boolean mask over ``df`` (needs observation_id and trait) of rows hit by a rule."""
    mask = pd.Series(False, index=df.index)
    r = rules[rules["scope"].isin(scopes)]
    if r.empty or df.empty:
        return mask
    c = df[["observation_id"]].merge(ctx[["observation_id", "dataset_id", "year"]], on="observation_id",
                                     how="left").set_index(df.index)
    for e in r.itertuples():
        hit = c["dataset_id"] == e.dataset_id
        if pd.notna(e.year):
            hit &= c["year"] == int(e.year)
        if pd.notna(e.trait):
            hit &= df["trait"] == e.trait
        mask |= hit.fillna(False)
    return mask


def _metrics(est: pd.Series, lo: pd.Series, hi: pd.Series, meas: pd.Series, relative: bool) -> dict:
    ok = np.isfinite(est) & np.isfinite(meas) & (meas > 0) & (est > 0)
    e, m = est[ok], meas[ok]
    n = int(ok.sum())
    if n < MIN_N:
        return {"n": n}
    le, lm = np.log(e), np.log(m)
    r = float(np.corrcoef(le, lm)[0, 1])
    rho = float(pd.Series(e).rank().corr(pd.Series(m).rank()))
    out = {"n": n, "r2_log": r**2, "spearman": rho}
    if not relative:
        slope = float(np.polyfit(lm, le, 1)[0])
        out.update({"rmse": float(np.sqrt(np.mean((e - m) ** 2))),
                    "nrmse": float(np.sqrt(np.mean((e - m) ** 2)) / np.mean(m)),
                    "bias_log": float(np.mean(le - lm)), "slope_log": slope,
                    "coverage90": float(np.mean((lo[ok] <= m) & (m <= hi[ok]))),
                    "ccc": float(_ccc(le, lm))})
    return out


def _ccc(x, y):
    sxy = np.cov(x, y)[0, 1]
    return 2 * sxy / (np.var(x, ddof=1) + np.var(y, ddof=1) + (np.mean(x) - np.mean(y)) ** 2)


# ---------------------------------------------------------------------------
# Part A
# ---------------------------------------------------------------------------


def validation_leaf(est: pd.DataFrame, td: pd.DataFrame, meas: pd.DataFrame, ctx: pd.DataFrame,
                    runs: pd.DataFrame):
    j = est.merge(meas, on=["observation_id", "trait"]).merge(
        ctx[["observation_id", "dataset_id", "sample_state", "foreoptic", "species_id",
             "growth_form"]], on="observation_id", how="left")
    # CRPS on the log scale (comparable across traits)
    tdm = td.merge(meas, on=["observation_id", "trait"])
    tdm = tdm[(tdm["value"] > 0) & (tdm["measured"] > 0)]
    crps = tdm.groupby(["observation_id", "run_id", "trait"]).apply(
        lambda g: crps_samples(np.log(g["value"].to_numpy()), np.log(g["measured"].iat[0])),
        include_groups=False).rename("crps_log").reset_index()
    j = j.merge(crps, on=["observation_id", "run_id", "trait"], how="left")
    j["covered90"] = (j["q05"] <= j["measured"]) & (j["measured"] <= j["q95"])
    rows = []
    for keys, g in j.groupby(["run_id", "trait", "sample_state"]):
        rows.append({"run_id": keys[0], "trait": keys[1], "sample_state": keys[2], "group": "all",
                     **_metrics(g["q50"], g["q05"], g["q95"], g["measured"], g["relative"].iat[0]),
                     "crps_log": g["crps_log"].mean()})
        for ds, h in g.groupby("dataset_id"):
            rows.append({"run_id": keys[0], "trait": keys[1], "sample_state": keys[2],
                         "group": f"dataset:{ds}",
                         **_metrics(h["q50"], h["q05"], h["q95"], h["measured"],
                                    h["relative"].iat[0]), "crps_log": h["crps_log"].mean()})
    metrics = pd.DataFrame(rows).pipe(with_runs, runs)
    return j.pipe(with_runs, runs), metrics


def validation_plant(td: pd.DataFrame, plant_meas: pd.DataFrame, ctx: pd.DataFrame) -> pd.DataFrame:
    """Plant-level (pooled-sample) traits vs the mean estimate over the plant's leaves."""
    if plant_meas.empty:
        return pd.DataFrame()
    t = td.merge(ctx[["observation_id", "plant_id"]], on="observation_id").dropna(subset=["plant_id"])
    pm = t.groupby(["plant_id", "run_id", "trait", "draw"])["value"].mean().reset_index()
    s = pm.groupby(["plant_id", "run_id", "trait"])["value"].agg(
        q50="median", q05=lambda x: x.quantile(0.05), q95=lambda x: x.quantile(0.95)).reset_index()
    return s.merge(plant_meas, on=["plant_id", "trait"])


def species_means(td: pd.DataFrame, meas: pd.DataFrame, ctx: pd.DataFrame, min_leaves=3):
    """Species-mean estimates with uncertainty propagated through the draws,
    vs species means of the measurements (same leaves)."""
    t = td.merge(ctx[["observation_id", "species_id", "dataset_id"]], on="observation_id")
    t = t.dropna(subset=["species_id"])
    m = meas.merge(ctx[["observation_id", "species_id"]], on="observation_id").dropna(subset=["species_id"])
    both = t.merge(m[["observation_id", "trait"]], on=["observation_id", "trait"])
    by_draw = both.groupby(["species_id", "run_id", "trait", "draw"])["value"].mean().reset_index()
    est = by_draw.groupby(["species_id", "run_id", "trait"])["value"].agg(
        q50="median", q05=lambda x: x.quantile(0.05), q95=lambda x: x.quantile(0.95)).reset_index()
    n = both.groupby(["species_id", "run_id", "trait"])["observation_id"].nunique().rename("n_leaves")
    est = est.merge(n.reset_index(), on=["species_id", "run_id", "trait"])
    mm = m.groupby(["species_id", "trait"])["measured"].agg(["mean", "std", "count"]).reset_index()
    mm.columns = ["species_id", "trait", "measured", "measured_sd", "n_measured"]
    out = est.merge(mm, on=["species_id", "trait"])
    return out[out["n_leaves"] >= min_leaves]


def spectral_residuals(spec: pd.DataFrame, ctx: pd.DataFrame, bin_nm=10.0):
    s = spec.merge(ctx[["observation_id", "sample_state", "foreoptic", "dataset_id"]],
                   on="observation_id", how="left")
    s["resid"] = s["y"] - s["yhat_mean"]
    s["covered90"] = (s["yhat_q05"] <= s["y"]) & (s["y"] <= s["yhat_q95"])
    s["wl_bin"] = (s["wavelength"] // bin_nm) * bin_nm
    byw = s.groupby(["run_id", "quantity", "fitted", "sample_state", "wl_bin"]).agg(
        resid_mean=("resid", "mean"), resid_sd=("resid", "std"),
        resid_q05=("resid", lambda x: x.quantile(0.05)), resid_q95=("resid", lambda x: x.quantile(0.95)),
        n=("resid", "size")).reset_index()
    per_leaf = s.groupby(["observation_id", "run_id", "quantity", "fitted"]).agg(
        rmse=("resid", lambda x: float(np.sqrt(np.mean(x**2)))), bias=("resid", "mean"),
        coverage90=("covered90", "mean")).reset_index()
    return byw, per_leaf


def heldout(per_leaf: pd.DataFrame, runs: pd.DataFrame) -> pd.DataFrame:
    """Out-of-sample prediction of quantities not used in the fit (e.g. T from R-only)."""
    h = per_leaf[~per_leaf["fitted"]].pipe(with_runs, runs)
    return h.groupby(["run_id", "model", "errmodel", "rt", "data_use", "quantity"]).agg(
        n=("observation_id", "nunique"), rmse=("rmse", "median"), bias=("bias", "median"),
        coverage90=("coverage90", "mean")).reset_index()


def intercompare(est: pd.DataFrame, runs: pd.DataFrame) -> pd.DataFrame:
    """Pairwise agreement between runs (different RTMs / error models) for the same trait."""
    e = est[~est["relative"]]
    wide = e.pivot_table(index=["observation_id", "trait"], columns="run_id", values="q50")
    rows = []
    for trait, w in wide.groupby(level="trait"):
        for a, b in combinations([c for c in w.columns if w[c].notna().any()], 2):
            x, y = w[a], w[b]
            ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
            if ok.sum() < MIN_N:
                continue
            lx, ly = np.log(x[ok]), np.log(y[ok])
            rows.append({"trait": trait, "run_a": a, "run_b": b, "n": int(ok.sum()),
                         "r2_log": float(np.corrcoef(lx, ly)[0, 1] ** 2),
                         "ccc_log": float(_ccc(lx, ly)),
                         "median_ratio": float(np.exp(np.median(ly - lx)))})
    return pd.DataFrame(rows)


def diagnostics(info: pd.DataFrame, summ: pd.DataFrame) -> pd.DataFrame:
    i = info[info["status"] == "ok"].copy()
    i["converged"] = (i["max_rhat"] < 1.01) & (i["min_ess_bulk"] >= 200) & (i["min_ess_tail"] >= 200)
    return i.groupby(["run_id", "model", "errmodel", "rt", "data_use"]).agg(
        n_fits=("observation_id", "size"), converged=("converged", "mean"),
        median_seconds=("seconds", "median"), p90_seconds=("seconds", lambda x: x.quantile(0.9)),
        total_hours=("seconds", lambda x: x.sum() / 3600),
        median_min_ess=("min_ess_bulk", "median"), max_rhat=("max_rhat", "max"),
        divergent_fits=("divergences", lambda x: float(np.mean(x > 0)))).reset_index()


USABLE = {"max_rhat": 1.05, "min_ess_bulk": 100, "divergent_share": 0.10}


def fit_screen(info: pd.DataFrame) -> pd.DataFrame:
    """Per fit: is it usable for trait analyses (R-hat <= 1.05, ESS_bulk >= 100,
    <= 10 % divergent transitions)? A looser screen than the ``converged`` rate in
    ``diagnostics``, which uses 1.01 / 200 for reporting. Failed fits are unusable."""
    i = info.copy()
    n = i["n_draws"] if "n_draws" in i else np.nan
    i["divergent_share"] = i["divergences"] / n
    ok = i["status"] == "ok"
    reasons = pd.Series("", index=i.index)
    reasons[~ok] = "error"
    for col, lim, label in (("max_rhat", USABLE["max_rhat"], "rhat"),
                            ("min_ess_bulk", USABLE["min_ess_bulk"], "ess"),
                            ("divergent_share", USABLE["divergent_share"], "divergences")):
        bad = ok & ((i[col] < lim) if col == "min_ess_bulk" else (i[col] > lim))
        reasons[bad] = reasons[bad].str.cat([label] * int(bad.sum()), sep=";").str.strip(";")
    i["usable"] = reasons == ""
    i["reason"] = reasons
    cols = ["observation_id", "run_id", "usable", "reason", "max_rhat", "min_ess_bulk", "divergent_share"]
    cols += [c for c in ("wall_margin", "wall_margin_q05", "frac_at_bound") if c in i]
    return i[cols]


def identifiability(summ: pd.DataFrame, priors: dict) -> pd.DataFrame:
    """Posterior / prior SD (log scale for positive parameters): < 1 = learned from data."""
    s = summ[summ["is_rtm"]].copy()

    def prior_sd_log(p):
        spec = priors["rtm"].get(p, {})
        if spec.get("dist") == "lognormal":
            return spec["sdlog"]
        if spec.get("dist") == "truncnormal":
            return spec["sigma"] / max(spec["mu"], 1e-9)
        if spec.get("dist") == "uniform":
            return (np.log(spec["upper"]) - np.log(spec["lower"])) / np.sqrt(12)
        return np.nan

    s["prior_sd_log"] = s["param"].map(prior_sd_log)
    s["post_sd_log"] = (np.log(s["q95"]) - np.log(s["q05"])) / (2 * 1.645)
    s["contraction"] = s["post_sd_log"] / s["prior_sd_log"]
    return s.groupby(["run_id", "model", "errmodel", "param"]).agg(
        contraction_median=("contraction", "median"),
        contraction_q90=("contraction", lambda x: x.quantile(0.9))).reset_index()


def posterior_correlations(draws: pd.DataFrame) -> pd.DataFrame:
    rows = []
    params = [c for c in draws.columns if c not in ("observation_id", "run_id", "draw", "model",
                                                     "errmodel", "rt", "data_use", "stage",
                                                     "only_rt", "sampler")]
    for (run_id, model), d in draws.groupby(["run_id", "model"]):
        ps = [p for p in params if d[p].notna().all()]
        cs = []
        for _, g in d.groupby("observation_id"):
            c = np.log(g[ps].clip(lower=1e-12)).corr()
            cs.append(c.to_numpy())
        med = np.nanmedian(np.stack(cs), axis=0)
        for i, a in enumerate(ps):
            for j, b in enumerate(ps):
                if i < j:
                    rows.append({"run_id": run_id, "model": model, "param_a": a, "param_b": b,
                                 "median_corr": med[i, j]})
    return pd.DataFrame(rows)


def error_model_comparison(metrics, held, diag, runs) -> pd.DataFrame:
    """One row per run: trait accuracy/calibration, held-out T, and cost."""
    m = metrics[(metrics["group"] == "all") & (metrics["sample_state"] == "fresh")]
    agg = m.groupby("run_id").agg(traits=("trait", "nunique"), coverage90=("coverage90", "mean"),
                                  crps_log=("crps_log", "mean"), r2_log=("r2_log", "mean"),
                                  nrmse=("nrmse", "mean")).reset_index()
    t = held[held["quantity"] == "T"][["run_id", "rmse", "coverage90"]].rename(
        columns={"rmse": "heldout_T_rmse", "coverage90": "heldout_T_coverage90"})
    out = runs.merge(agg, on="run_id", how="left").merge(t, on="run_id", how="left").merge(
        diag[["run_id", "converged", "median_seconds", "total_hours"]], on="run_id", how="left")
    return out


def data_use_comparison(leaf_val: pd.DataFrame) -> pd.DataFrame:
    """R-only vs all-quantity fits of the same model/error model on the same leaves."""
    rows = []
    ronly = leaf_val[leaf_val["data_use"] == "R"]
    for (model, err), g in ronly.groupby(["model", "errmodel"]):
        full = leaf_val[(leaf_val["model"] == model) & (leaf_val["errmodel"] == err)
                        & (leaf_val["data_use"] == "all") & ~leaf_val["rt"].astype(bool)]
        j = g.merge(full, on=["observation_id", "trait"], suffixes=("_R", "_RT"))
        for t, h in j.groupby("trait"):
            if len(h) < MIN_N or h["relative_R"].iat[0]:
                continue
            row = {"model": model, "errmodel": err, "trait": t, "n": len(h)}
            for suf in ("_R", "_RT"):
                m = _metrics(h["q50" + suf], h["q05" + suf], h["q95" + suf], h["measured" + suf], False)
                row.update({f"{k}{suf}": v for k, v in m.items() if k != "n"})
                row[f"crps_log{suf}"] = h["crps_log" + suf].mean()
            rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Part B
# ---------------------------------------------------------------------------


def best_config(leaf_val: pd.DataFrame) -> pd.DataFrame:
    """Per trait x foreoptic class: the run with the lowest mean log-CRPS on
    measured fresh leaves, with leave-one-dataset-out (LODO) evaluation of the
    selection procedure itself."""
    v = leaf_val[(leaf_val["sample_state"] == "fresh") & ~leaf_val["relative"]].dropna(subset=["crps_log"])
    v = v.assign(optic=v["foreoptic"].map({"integrating_sphere": "sphere", "leaf_clip": "clip",
                                            "contact_probe": "clip"}).fillna("other"))
    rows = []
    for (trait, optic), g in list(v.groupby(["trait", "optic"])) + [((t, "all"), g) for t, g in v.groupby("trait")]:
        # only runs fitted to every leaf in the group compete (fair comparison)
        cover = g.groupby("run_id")["observation_id"].nunique()
        full = cover[cover == cover.max()].index
        g = g[g["run_id"].isin(full)]
        score = g.groupby("run_id")["crps_log"].mean().sort_values()
        if score.empty:
            continue
        lodo = []
        for ds in g["dataset_id"].unique():
            train, test = g[g["dataset_id"] != ds], g[g["dataset_id"] == ds]
            if train.empty:
                continue
            pick = train.groupby("run_id")["crps_log"].mean().idxmin()
            lodo.append(test.loc[test["run_id"] == pick, "crps_log"].mean())
        rows.append({"trait": trait, "optic": optic, "best_run": score.index[0],
                     "crps_log": score.iloc[0], "runner_up": score.index[1] if len(score) > 1 else None,
                     "lodo_crps_log": float(np.nanmean(lodo)) if lodo else np.nan,
                     "n_leaves": int(g["observation_id"].nunique()),
                     "scores": json.dumps(score.round(4).to_dict())})
    return pd.DataFrame(rows)


def best_estimates(td: pd.DataFrame, best: pd.DataFrame, ctx: pd.DataFrame) -> pd.DataFrame:
    """Draws of the selected configuration per observation and trait."""
    c = ctx[["observation_id", "foreoptic"]].assign(
        optic=lambda d: d["foreoptic"].map({"integrating_sphere": "sphere", "leaf_clip": "clip",
                                            "contact_probe": "clip"}).fillna("other"))
    sel = best.set_index(["trait", "optic"])["best_run"].to_dict()
    sel_all = best[best["optic"] == "all"].set_index("trait")["best_run"].to_dict()
    t = td.merge(c[["observation_id", "optic"]], on="observation_id")
    choice = [sel.get((tr, o), sel_all.get(tr)) for tr, o in zip(t["trait"], t["optic"])]
    return t[t["run_id"].to_numpy() == np.array(choice, dtype=object)].drop(columns="optic")


def variance_partition(values: pd.DataFrame, ctx: pd.DataFrame, source: str) -> pd.DataFrame:
    """Variance of log trait among datasets, among species within datasets,
    and within species (method of moments, balanced-ish nested ANOVA)."""
    v = values.merge(ctx[["observation_id", "species_id", "dataset_id"]], on="observation_id")
    v = v.dropna(subset=["species_id"])
    v = v[v["value"] > 0]
    rows = []
    for trait, g in v.groupby("trait"):
        g = g.assign(y=np.log(g["value"]))
        sp_n = g.groupby("species_id")["y"].transform("size")
        g = g[sp_n >= 2]
        if g["species_id"].nunique() < 3:
            continue
        total = g["y"].var()
        within = g.groupby("species_id")["y"].var().mean()  # mean within-species variance
        sp_means = g.groupby(["dataset_id", "species_id"])["y"].mean()
        among_ds = sp_means.groupby(level=0).mean().var() if sp_means.index.get_level_values(0).nunique() > 1 else 0.0
        among_sp = max(total - within - (among_ds if np.isfinite(among_ds) else 0), 0)
        rows.append({"trait": trait, "source": source, "n_leaves": len(g),
                     "n_species": g["species_id"].nunique(), "total_var": total,
                     "within_species": within / total, "among_species": among_sp / total,
                     "among_datasets": (among_ds if np.isfinite(among_ds) else 0) / total})
    return pd.DataFrame(rows)


def treatment_effects(values: pd.DataFrame, db: dict, source: str) -> pd.DataFrame:
    """Log response ratio of each treatment level vs control, per dataset x factor,
    with bootstrap 90 % intervals."""
    links = db["observation_treatments"].merge(db["treatments"], on="treatment_id")
    links = links[links["is_control"].notna()]
    rng = np.random.default_rng(1)
    rows = []
    v = values.merge(links[["observation_id", "dataset_id", "factor", "level", "is_control"]],
                     on="observation_id")
    v = v[v["value"] > 0]
    for (ds, factor, trait), g in v.groupby(["dataset_id", "factor", "trait"]):
        ctl = np.log(g.loc[g["is_control"].astype(bool), "value"].to_numpy())
        if len(ctl) < 3:
            continue
        for level, h in g[~g["is_control"].astype(bool)].groupby("level"):
            x = np.log(h["value"].to_numpy())
            if len(x) < 3:
                continue
            boot = [rng.choice(x, len(x)).mean() - rng.choice(ctl, len(ctl)).mean() for _ in range(500)]
            rows.append({"dataset_id": ds, "factor": factor, "level": level, "trait": trait,
                         "source": source, "lrr": x.mean() - ctl.mean(),
                         "lrr_q05": np.quantile(boot, 0.05), "lrr_q95": np.quantile(boot, 0.95),
                         "n_treat": len(x), "n_control": len(ctl)})
    return pd.DataFrame(rows)


def correlations(best_sum: pd.DataFrame, meas: pd.DataFrame, ctx: pd.DataFrame, min_species_n=10):
    """Leaf-level and species-mean Spearman correlations among RTM estimates
    and measured traits, plus within-species correlations."""
    e = best_sum.assign(var="rtm:" + best_sum["trait"])[["observation_id", "var", "q50"]].rename(
        columns={"q50": "value"})
    m = meas.assign(var="meas:" + meas["trait"])[["observation_id", "var", "measured"]].rename(
        columns={"measured": "value"})
    long = pd.concat([e, m]).merge(ctx[["observation_id", "species_id"]], on="observation_id")
    wide = long.pivot_table(index=["observation_id", "species_id"], columns="var", values="value")
    keep = [c for c in wide.columns if wide[c].notna().sum() >= 20]
    wide = np.log(wide[keep].where(wide[keep] > 0))
    rows, within = [], []
    sp_means = wide.groupby(level="species_id").mean()
    sp_counts = wide.groupby(level="species_id").count()
    for a, b in combinations(keep, 2):
        ok = wide[[a, b]].dropna()
        if len(ok) >= 20:
            rows.append({"var_a": a, "var_b": b, "level": "leaf", "n": len(ok),
                         "spearman": ok[a].corr(ok[b], method="spearman")})
        okm = sp_means[[a, b]][(sp_counts[a] >= 3) & (sp_counts[b] >= 3)].dropna()
        if len(okm) >= 8:
            rows.append({"var_a": a, "var_b": b, "level": "species_mean", "n": len(okm),
                         "spearman": okm[a].corr(okm[b], method="spearman")})
        for sid, g in ok.groupby(level="species_id"):
            if len(g) >= min_species_n:
                within.append({"var_a": a, "var_b": b, "species_id": sid, "n": len(g),
                               "spearman": g[a].corr(g[b], method="spearman")})
    return pd.DataFrame(rows), pd.DataFrame(within)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(argv=None):
    from leafrtm.run import load_config

    ap = argparse.ArgumentParser()
    ap.add_argument("--fitset-dir", required=True)
    a = ap.parse_args(argv)
    acfg = load_config().get("analysis", {}) or {}
    d = Path(a.fitset_dir)
    out = d / "analysis"
    out.mkdir(parents=True, exist_ok=True)
    S = {p: pd.read_parquet(d / "summary" / f"{p}.parquet") for p in ("summary", "draws", "spectra", "info")}
    runs = S["info"][["run_id", "model", "errmodel", "rt", "data_use"]].drop_duplicates("run_id")
    db = load()
    ctx = observation_context(db)
    screen = fit_screen(S["info"])
    usable = screen.loc[screen["usable"], ["observation_id", "run_id"]]
    # trait analyses use usable fits only; spectral fit tables keep every fit
    td = trait_draws(S["draws"].merge(usable, on=["observation_id", "run_id"]))
    est = summarize_draws(td)
    meas_all = measured_traits(db)
    tables = {}
    # suspect groups are found on every measurement; exclusions with scope "all" then
    # remove measurements from validation and best-configuration selection too
    leaf_all, _ = validation_leaf(est, td, meas_all, ctx, runs)
    suspects = suspect_groups(leaf_all, ctx, acfg.get("flag"))
    rules = exclusion_rules(acfg, suspects)
    tables["suspect_groups"], tables["exclusion_rules"] = suspects, rules
    meas = meas_all[~excluded(meas_all, ctx, rules, scopes=("all",))]
    leaf, metrics = validation_leaf(est, td, meas, ctx, runs)
    tables["estimates"] = est.pipe(with_runs, runs)
    tables["validation_leaf"] = leaf
    tables["validation_metrics"] = metrics
    tables["validation_plant"] = validation_plant(td, plant_traits(db), ctx)
    tables["validation_species"] = species_means(td, meas, ctx).pipe(with_runs, runs)
    byw, per_leaf = spectral_residuals(S["spectra"], ctx)
    tables["spectral_residuals"] = byw.pipe(with_runs, runs)
    tables["spectral_fit_leaf"] = per_leaf.pipe(with_runs, runs)
    tables["heldout"] = heldout(per_leaf, runs)
    tables["data_use_comparison"] = data_use_comparison(leaf)
    tables["intercompare"] = intercompare(est, runs)
    tables["diagnostics"] = diagnostics(S["info"], S["summary"]).merge(
        screen.groupby("run_id")["usable"].mean().rename("usable").reset_index(), on="run_id")
    tables["fit_screen"] = screen.pipe(with_runs, runs)
    tables["excluded_fits"] = tables["fit_screen"][~tables["fit_screen"]["usable"]]
    tables["identifiability"] = identifiability(S["summary"], load_priors())
    tables["posterior_correlations"] = posterior_correlations(S["draws"])
    tables["error_model_comparison"] = error_model_comparison(metrics, tables["heldout"],
                                                              tables["diagnostics"], runs)
    best = best_config(leaf)
    tables["best_config"] = best
    bd = best_estimates(td, best, ctx)
    bsum = summarize_draws(bd) if len(bd) else pd.DataFrame(columns=["observation_id", "trait", "q50"])
    tables["best_estimates"] = bsum
    # ecological analyses: drop every excluded group (estimates and measurements alike)
    bsum = bsum[~excluded(bsum, ctx, rules)]
    vals_rtm = bsum[["observation_id", "trait", "q50"]].rename(columns={"q50": "value"})
    meas_eco = meas[~excluded(meas, ctx, rules)]
    vals_meas = meas_eco.rename(columns={"measured": "value"})
    tables["exclusions_applied"] = (
        pd.concat([tables["best_estimates"].assign(kind="rtm"), meas_all.assign(kind="measured")])
        .loc[lambda d: excluded(d, ctx, rules)]
        .merge(ctx[["observation_id", "dataset_id", "year"]], on="observation_id")
        .groupby(["dataset_id", "year", "trait", "kind"]).size().rename("n").reset_index())
    same = vals_rtm.merge(vals_meas[["observation_id", "trait"]], on=["observation_id", "trait"])
    tables["variance"] = pd.concat([
        variance_partition(vals_rtm, ctx, "rtm_all_leaves"),
        variance_partition(same, ctx, "rtm_measured_leaves"),
        variance_partition(vals_meas, ctx, "measured")], ignore_index=True)
    tables["treatments"] = pd.concat([treatment_effects(vals_rtm, db, "rtm"),
                                      treatment_effects(vals_meas, db, "measured")], ignore_index=True)
    corr, within = correlations(bsum, meas_eco, ctx)
    tables["correlations"], tables["within_species_correlations"] = corr, within
    for name, df in tables.items():
        df.to_parquet(out / f"{name}.parquet", index=False)
        print(f"{name}: {len(df)}")


if __name__ == "__main__":
    main()
