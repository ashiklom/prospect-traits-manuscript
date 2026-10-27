"""Trait estimates from posterior draws, joined with measured traits."""

from __future__ import annotations

from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]


@cache
def mapping() -> dict:
    with open(ROOT / "config" / "rtm_traits.yaml") as f:
        return yaml.safe_load(f)


def trait_draws(draws: pd.DataFrame) -> pd.DataFrame:
    """Long table: observation_id, run_id, model, trait, draw, value."""
    out = []
    for (run_id, model), d in draws.groupby(["run_id", "model"]):
        for m in mapping().get(model, []):
            expr = m.get("expr", m.get("param"))
            try:
                v = d.eval(expr, engine="python") * m["factor"]
            except Exception:
                continue
            out.append(pd.DataFrame({
                "observation_id": d["observation_id"].to_numpy(), "run_id": run_id,
                "model": model, "trait": m["trait"], "draw": d["draw"].to_numpy(),
                "value": np.asarray(v, float), "relative": bool(m.get("relative", False))}))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def summarize_draws(td: pd.DataFrame) -> pd.DataFrame:
    g = td.groupby(["observation_id", "run_id", "model", "trait", "relative"])["value"]
    s = g.agg(mean="mean", sd="std", q05=lambda x: x.quantile(0.05), q50="median",
              q95=lambda x: x.quantile(0.95)).reset_index()
    return s


def measured_traits(db: dict) -> pd.DataFrame:
    """Leaf-level measured traits per observation (plus lignin+cellulose per area)."""
    tr = db["traits"]
    leaf = tr[tr["trait_level"] == "leaf"].dropna(subset=["observation_id"])
    wide = leaf.pivot_table(index="observation_id", columns="trait", values="value", aggfunc="mean")
    if {"lignin_mass", "cellulose_mass", "lma"} <= set(wide):
        wide["lignin_cellulose_area"] = (wide["lignin_mass"] + wide["cellulose_mass"]) / 100 * wide["lma"]
    return wide.reset_index().melt(id_vars="observation_id", var_name="trait",
                                   value_name="measured").dropna()


def plant_traits(db: dict) -> pd.DataFrame:
    """Traits measured on pooled samples (plant level), per plant_id."""
    tr = db["traits"]
    pl = tr[tr["trait_level"] == "plant"].merge(db["entities"][["entity_id", "plant_id"]],
                                               on="entity_id")
    return pl.groupby(["plant_id", "trait"], as_index=False)["value"].mean().rename(
        columns={"value": "measured"})


def crps_samples(x: np.ndarray, y: float) -> float:
    """Empirical CRPS of a sample ``x`` for observation ``y``."""
    x = np.sort(np.asarray(x, float))
    n = len(x)
    term1 = np.mean(np.abs(x - y))
    # E|X - X'| via sorted samples
    term2 = 2 * np.sum(x * (2 * np.arange(1, n + 1) - n - 1)) / n**2
    return term1 - 0.5 * term2
