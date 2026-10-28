"""End-to-end inversion checks (short chains)."""

import numpy as np
import pytest

from leafrtm import likelihood as L
from leafrtm.fit import FitSpec, Inverter, data_arrays
from leafrtm.obsop import Block, make_obsdata
from leafrtm.rtm import get_model

TRUTH = dict(N=1.6, Cab=45, Car=9, Canth=1.5, Cbrown=0.05, Cw=0.012, Cm=0.006)


def _leaf(step, seed, quantities=("R", "T")):
    rng = np.random.default_rng(seed)
    m = get_model("prospectd")
    R, T = m.simulate(TRUTH)
    wl = np.arange(400, 2501, step)
    i = (wl - 400).astype(int)
    blocks = [Block(q, wl, {"R": R, "T": T}[q][i] + rng.normal(0, 0.004, len(wl)))
              for q in quantities]
    return make_obsdata(blocks, m.wavelength)


@pytest.fixture(scope="module")
def ou_inverter():
    spec = FitSpec("prospectd", "ou", True, ("R", "T"))
    return spec, Inverter(spec, data_arrays(_leaf(1.0, 0), spec))


def test_swapped_data_equals_fresh_compile(ou_inverter):
    spec, inv = ou_inverter
    B = data_arrays(_leaf(3.0, 1), spec)
    tr1, _ = inv.sample(B, draws=100, tune=150, chains=1, seed=5)
    tr2, _ = Inverter(spec, B).sample(B, draws=100, tune=150, chains=1, seed=5)
    np.testing.assert_array_equal(tr1.posterior["Cab"].values, tr2.posterior["Cab"].values)


def test_recovers_truth(ou_inverter):
    spec, inv = ou_inverter
    tr, _ = inv.sample(data_arrays(_leaf(2.0, 2), spec), draws=300, tune=400, chains=2, seed=3)
    for p in ("N", "Cab", "Cw", "Cm"):
        x = tr.posterior[p].values.ravel()
        lo, hi = np.quantile(x, [0.001, 0.999])
        assert lo * 0.97 <= TRUTH[p] <= hi * 1.03, (p, lo, hi)


def test_ss_model_runs():
    spec = FitSpec("leafsip", "ss", False, ("R",))
    od = _leaf(5.0, 3, quantities=("R",))
    nb, b, S = L.vecchia(od.wavelength, od.start, lambda d: L.matern32(d, 40.0, 0.02))
    arrays = data_arrays(od, spec, {"nb": nb, "b": b, "S": S})
    tr, _ = Inverter(spec, arrays).sample(arrays, draws=50, tune=100, chains=1, seed=1)
    assert np.isfinite(tr.posterior["Cab"].values).all()


def test_error_priors_sample_without_divergences():
    """Every error-model prior must map the whole unconstrained space into its support:
    a truncated lognormal built with lower=None admitted ell <= 0 (17 % divergences)."""
    import nutpie
    import pymc as pm

    from leafrtm.fit import _prior, load_priors

    with pm.Model() as model:
        for name, spec in load_priors()["error"].items():
            _prior(name, spec, shape=2)
    trace = nutpie.sample(nutpie.compile_pymc_model(model, backend="numba"), draws=500, tune=500,
                          chains=2, seed=1, progress_bar=False)
    assert int(np.asarray(trace.sample_stats["diverging"]).sum()) == 0
