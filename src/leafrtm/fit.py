"""Bayesian inversion of one leaf observation with PyMC + nutpie.

A fit is specified by an RTM, an error model (``ou``, ``ss``, ``ld``), whether
R/T errors are correlated (``rt``), and the set of observed quantities. The
PyMC model is built with ``pm.Data`` containers, compiled once per such
specification (``Inverter``), and re-used across leaves with
``CompiledModel.with_data``; compilation is ~10 s, while a fit is seconds to
minutes.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path

import numpy as np
import pymc as pm
import pytensor.tensor as pt
import yaml

from leafrtm import likelihood as L
from leafrtm.obsop import QUANTITIES, ObsData
from leafrtm.rtm import get_model

CONFIG = Path(__file__).resolve().parents[2] / "config"
ERRMODELS = ("ou", "ss", "ld")


@cache
def load_priors(path: str | None = None) -> dict:
    with open(path or CONFIG / "priors.yaml") as f:
        return yaml.safe_load(f)


def _prior(name, spec, **kw):
    d = spec["dist"]
    if d == "lognormal":
        return pm.LogNormal(name, mu=spec["meanlog"], sigma=spec["sdlog"], **kw)
    if d == "truncnormal":
        return pm.TruncatedNormal(name, mu=spec["mu"], sigma=spec["sigma"],
                                  lower=spec.get("lower"), upper=spec.get("upper"), **kw)
    if d == "uniform":
        return pm.Uniform(name, lower=spec["lower"], upper=spec["upper"], **kw)
    if d == "halfnormal":
        return pm.HalfNormal(name, sigma=spec["sigma"], **kw)
    if d == "beta":
        return pm.Beta(name, alpha=spec["alpha"], beta=spec["beta"], **kw)
    raise ValueError(f"Unknown prior {d}")


def prior_median(spec) -> float:
    d = spec["dist"]
    if d == "lognormal":
        return float(np.exp(spec["meanlog"]))
    if d == "truncnormal":
        return float(np.clip(spec["mu"], spec.get("lower", -np.inf), spec.get("upper", np.inf)))
    if d == "uniform":
        return 0.5 * (spec["lower"] + spec["upper"])
    raise ValueError(d)


@dataclass(frozen=True)
class FitSpec:
    model: str
    errmodel: str
    rt: bool
    quantities: tuple[str, ...]

    def __post_init__(self):
        if self.errmodel not in ERRMODELS:
            raise ValueError(self.errmodel)
        if self.rt and not {"R", "T"} <= set(self.quantities):
            object.__setattr__(self, "rt", False)

    @property
    def key(self) -> str:
        return f"{self.model}-{self.errmodel}{'-rt' if self.rt else ''}-{''.join(self.quantities)}"


def data_arrays(od: ObsData, spec: FitSpec, extra: dict | None = None) -> dict:
    """Flat arrays fed to ``pm.Data`` (shapes vary by leaf)."""
    local = {QUANTITIES.index(q): i for i, q in enumerate(spec.quantities)}
    if set(od.quantities) != set(spec.quantities):
        raise ValueError(f"Observation has {od.quantities}, spec expects {spec.quantities}")
    lq = np.array([local[q] for q in od.qcode], dtype=np.int64)
    out = {"y": od.y, "flat_idx": od.flat_idx, "w": od.w, "lq": lq,
           "dlam": od.dlam, "start": od.start.astype(np.int8), "pair": od.pair}
    if spec.errmodel in ("ss", "ld"):
        missing = {"nb", "b", "S"} - set(extra or {})
        if missing:
            raise ValueError(f"{spec.errmodel} needs precomputed {sorted(missing)}")
    out.update(extra or {})
    return out


def build_model(spec: FitSpec, arrays: dict, priors: dict | None = None) -> pm.Model:
    priors = priors or load_priors()
    rtm = get_model(spec.model)
    nq = len(spec.quantities)
    with pm.Model() as model:
        D = {k: pm.Data(k, v) for k, v in arrays.items()}
        theta = {p: _prior(p, priors["rtm"][p]) for p in rtm.params}
        R, T = rtm.simulate_pt(theta)
        with np.errstate(divide="ignore"):
            Q = pt.stack([R, T, 1 - R - T, -pt.log10(R)]).ravel()
        yhat = pt.sum(D["w"] * Q[D["flat_idx"]], axis=1)
        pm.Deterministic("max_pred", pt.max(yhat))  # cheap sanity trace
        e = D["y"] - yhat
        d = {"qcode": D["lq"], "dlam": D["dlam"], "start": pt.neq(D["start"], 0), "pair": D["pair"]}
        rho = None
        if spec.rt and spec.errmodel != "ld":
            rho = pm.Deterministic("rho", 2 * _prior("rho_raw", priors["error"]["rho"]) - 1)
        if spec.errmodel == "ou":
            s0 = _prior("sigma0", priors["error"]["sigma0"], shape=nq)
            s1 = _prior("sigma1", priors["error"]["sigma1"], shape=nq)
            ell = _prior("ell", priors["error"]["ell"], shape=nq)
            ll = L.loglik_ou(e, yhat, d, s0, s1, ell, rho=rho)
        elif spec.errmodel == "ss":
            s0 = _prior("sigma0", priors["error"]["sigma0"], shape=nq)
            s1 = _prior("sigma1", priors["error"]["sigma1"], shape=nq)
            d.update(nb=D["nb"], b=D["b"], S=D["S"])
            ll = L.loglik_ss(e, yhat, d, s0, s1, rho=rho)
        else:
            kappa = _prior("kappa", priors["error"]["kappa"], shape=nq)
            d.update({k: D[k] for k in ("mu", "sqrtD", "nb", "b", "S", "W", "Minv", "logdetM")})
            ll = L.loglik_ld(e, d, kappa)
        pm.Potential("loglik", ll)
    return model


def initial_point(spec: FitSpec, priors: dict | None = None) -> dict:
    priors = priors or load_priors()
    rtm = get_model(spec.model)
    return {p: prior_median(priors["rtm"][p]) for p in rtm.params}


@dataclass
class Inverter:
    """Compiled sampler for one FitSpec, re-used across leaves."""

    spec: FitSpec
    template: dict
    priors: dict | None = None
    compiled: object = field(default=None, repr=False)
    compile_time: float = 0.0

    def __post_init__(self):
        import nutpie

        t0 = time.perf_counter()
        model = build_model(self.spec, self.template, self.priors)
        self.compiled = nutpie.compile_pymc_model(
            model, backend="numba", initial_points=initial_point(self.spec, self.priors))
        self.compile_time = time.perf_counter() - t0

    def sample(self, arrays: dict, draws=1000, tune=1000, chains=4, seed=0, **kwargs):
        import nutpie

        used = {v.name for v in self.compiled.shared_var_keys}
        cm = self.compiled.with_data(**{k: v for k, v in arrays.items() if k in used})
        t0 = time.perf_counter()
        trace = nutpie.sample(cm, draws=draws, tune=tune, chains=chains, cores=1, seed=seed,
                              progress_bar=False, **kwargs)
        return trace, time.perf_counter() - t0
