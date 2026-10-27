"""Forward-model tests for all leaf RTMs."""

import numpy as np
import pytensor
import pytensor.tensor as pt
import pytest

from leafrtm.rtm import MODELS, get_model
from leafrtm.rtm.e1 import e1_fast_np
from scipy.special import exp1

# Representative parameter sets per model (typical, thick, thin/pale).
CASES = {
    "prospect4": [dict(N=1.5, Cab=40, Cw=0.01, Cm=0.009), dict(N=2.5, Cab=90, Cw=0.03, Cm=0.02),
                  dict(N=1.1, Cab=1, Cw=0.002, Cm=0.001)],
    "prospect5": [dict(N=1.5, Cab=40, Car=10, Cw=0.01, Cm=0.009),
                  dict(N=1.1, Cab=1, Car=0.1, Cw=0.002, Cm=0.001)],
    "prospect5b": [dict(N=1.5, Cab=40, Car=10, Cbrown=0.1, Cw=0.01, Cm=0.009),
                   dict(N=1.8, Cab=2, Car=6, Cbrown=0.8, Cw=0.006, Cm=0.004)],
    "prospectd": [dict(N=1.5, Cab=40, Car=10, Canth=2, Cbrown=0.1, Cw=0.01, Cm=0.009),
                  dict(N=2.5, Cab=90, Car=20, Canth=5, Cbrown=0.5, Cw=0.03, Cm=0.02),
                  dict(N=1.1, Cab=1, Car=0.1, Canth=0.01, Cbrown=0.001, Cw=0.002, Cm=0.001),
                  dict(N=1.5, Cab=400, Car=100, Canth=30, Cbrown=3, Cw=0.08, Cm=0.05)],
}


def test_every_model_has_cases():
    assert set(CASES) == set(MODELS)


@pytest.mark.parametrize("name", MODELS)
def test_physical_range(name):
    m = get_model(name)
    for case in CASES[name]:
        R, T = m.simulate(case)
        assert R.shape == T.shape == m.wavelength.shape
        assert np.all(np.isfinite(R)) and np.all(np.isfinite(T))
        assert np.all((R >= 0) & (T >= 0) & (R + T <= 1 + 1e-12))


def test_fast_e1_accuracy():
    x = np.geomspace(1e-6, 700, 5000)
    np.testing.assert_allclose(e1_fast_np(x), exp1(x), rtol=1e-11)


def test_prospectd_reference_values():
    # prosail (Gomez-Dans) PROSPECT-D values from rrtm/tests/testthat/test.prospect.R
    wl = np.array([450, 550, 680, 800, 1450, 1940, 2200]) - 400
    ref_r = [0.0411123436, 0.1115443677, 0.0359676508, 0.4351475168, 0.1650296676,
             0.0373651153, 0.1547468978]
    ref_t = [0.0007884939, 0.1060743066, 0.0051708164, 0.4670787360, 0.2096989879,
             0.0489542338, 0.2531362617]
    m = get_model("prospectd")
    for e1 in ("exact", "fast"):
        R, T = m.simulate(CASES["prospectd"][0], e1=e1)
        np.testing.assert_allclose(R[wl], ref_r, atol=1e-9)
        np.testing.assert_allclose(T[wl], ref_t, atol=1e-9)


def _prosail(version, case, km=None, cm=None):
    prosail = pytest.importorskip("prosail")
    _, R, T = prosail.run_prospect(case["N"], case["Cab"], case.get("Car", 0), case.get("Cbrown", 0),
                                   case["Cw"], case["Cm"] if cm is None else cm,
                                   ant=case.get("Canth", 0), prospect_version=version, km=km)
    return R, T


@pytest.mark.parametrize("case", CASES["prospectd"])
def test_prospectd_matches_prosail(case):
    R, T = get_model("prospectd").simulate(case, e1="exact")
    Rp, Tp = _prosail("D", case)
    np.testing.assert_allclose(R, Rp, atol=1e-8)
    np.testing.assert_allclose(T, Tp, atol=1e-8)


@pytest.mark.parametrize("case", CASES["prospect5b"])
def test_prospect5b_matches_prosail(case):
    # prosail stores its PROSPECT-5 table in float32 (coefficients agree to
    # ~6e-8 relative), hence the looser tolerance than for PROSPECT-D.
    R, T = get_model("prospect5b").simulate(case, e1="exact")
    Rp, Tp = _prosail("5", case)
    np.testing.assert_allclose(R, Rp, atol=1e-7)
    np.testing.assert_allclose(T, Tp, atol=1e-7)


def test_prospect5_is_5b_without_brown():
    c5 = CASES["prospect5"][0]
    np.testing.assert_allclose(get_model("prospect5").simulate(c5)[0],
                               get_model("prospect5b").simulate({**c5, "Cbrown": 0.0})[0], atol=1e-14)


def _compile(m):
    sym = {p: pt.scalar(p) for p in m.params}
    R, T = m.simulate_pt(sym)
    g = pt.grad(R.sum() + 0.5 * T.sum(), list(sym.values()))
    return pytensor.function(list(sym.values()), [R, T, *g])


@pytest.mark.parametrize("name", MODELS)
def test_pytensor_matches_numpy_and_gradients(name):
    m = get_model(name)
    f = _compile(m)
    for case in CASES[name]:
        args = [case[p] for p in m.params]
        R, T, *grads = f(*args)
        Rn, Tn = m.simulate(case)
        np.testing.assert_allclose(R, Rn, atol=1e-12)
        np.testing.assert_allclose(T, Tn, atol=1e-12)
        for p, g in zip(m.params, grads):
            h = 1e-6 * max(abs(case[p]), 1e-3)
            up = m.simulate({**case, p: case[p] + h})
            dn = m.simulate({**case, p: case[p] - h})
            fd = ((up[0] - dn[0]).sum() + 0.5 * (up[1] - dn[1]).sum()) / (2 * h)
            assert g == pytest.approx(fd, rel=1e-4, abs=1e-6), (name, p)
