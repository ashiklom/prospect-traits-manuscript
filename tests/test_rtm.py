"""Forward-model tests for all leaf RTMs."""

import numpy as np
import pytensor
import pytensor.tensor as pt
import pytest

from leafrtm.rtm import MODELS, get_model
from leafrtm.rtm.e1 import e1_fast_np
from leafrtm.rtm.liberty import ABSORBERS
from liberty_rrtm import liberty_rrtm
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
    "prospectpro": [dict(N=1.5, Cab=40, Car=10, Canth=2, Cbrown=0.1, Cw=0.01, Cprot=0.001, Ccbc=0.008),
                    dict(N=2.0, Cab=60, Car=12, Canth=1, Cbrown=0.0, Cw=0.02, Cprot=0.002, Ccbc=0.015)],
    "liberty": [dict(xu=0.045, thick=1.6, baseline=0.0006, element=2, c_factor=200, w_factor=100,
                     l_factor=40, p_factor=1),
                dict(xu=0.02, thick=3.4, baseline=0.0004, element=1, c_factor=400, w_factor=250,
                     l_factor=60, p_factor=1.5)],
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


@pytest.mark.parametrize("case", CASES["prospectpro"])
def test_prospectpro_matches_prosail_with_combined_dry_matter(case):
    # PRO = PROSPECT-D with the dry-matter term split into protein + CBC; an
    # independent check passes K_prot*Cprot + K_cbc*Ccbc to prosail as km with cm = 1.
    m = get_model("prospectpro")
    R, T = m.simulate(case, e1="exact")
    K = dict(zip(m.constituents, m.K.T))
    km = K["Cprot"] * case["Cprot"] + K["Ccbc"] * case["Ccbc"]
    Rp, Tp = _prosail("D", {**case, "Cm": 1.0}, km=km)
    np.testing.assert_allclose(R, Rp, atol=1e-8)
    np.testing.assert_allclose(T, Tp, atol=1e-8)


def test_prospectpro_argument_order():
    # Regression for the old Julia code, which swapped Cprot and Ccbc.
    m = get_model("prospectpro")
    case = CASES["prospectpro"][0]
    swapped = {**case, "Cprot": case["Ccbc"], "Ccbc": case["Cprot"]}
    assert not np.allclose(m.simulate(case)[0], m.simulate(swapped)[0])
    K = dict(zip(m.constituents, m.K.T))
    assert K["Cprot"][0] == 0  # protein absorbs only in the NIR/SWIR


def test_prospect5_is_5b_without_brown():
    c5 = CASES["prospect5"][0]
    np.testing.assert_allclose(get_model("prospect5").simulate(c5)[0],
                               get_model("prospect5b").simulate({**c5, "Cbrown": 0.0})[0], atol=1e-14)


@pytest.mark.parametrize("case", CASES["liberty"])
def test_liberty_matches_converged_rrtm(case):
    m = get_model("liberty")
    R, T, Rinf = m.simulate(case, return_rinf=True)
    Ro, To, Rinfo = liberty_rrtm(m.K, m.fixed["D"], **case, n_iter=20000)
    np.testing.assert_allclose(R, Ro, atol=1e-11)
    np.testing.assert_allclose(T, To, atol=1e-11)
    np.testing.assert_allclose(Rinf, Rinfo, atol=1e-11)
    # rrtm's 50 iterations are close to, but not exactly at, the fixed point
    Ro50, To50, _ = liberty_rrtm(m.K, m.fixed["D"], **case)
    assert np.abs(R - Ro50).max() < 1e-4


@pytest.mark.parametrize("thick", [1.2, 2.0, 2.7, 3.0, 4.0, 6.5])
def test_liberty_integer_layers_match_loop(thick):
    m = get_model("liberty")
    case = {**CASES["liberty"][0], "thick": thick}
    R, T = m.simulate(case)
    Ro, To, _ = liberty_rrtm(m.K, m.fixed["D"], **case, n_iter=20000)
    np.testing.assert_allclose(R, Ro, atol=1e-11)
    np.testing.assert_allclose(T, To, atol=1e-11)


@pytest.mark.parametrize("k", [2, 3, 4])
def test_liberty_continuous_across_integer_thickness(k):
    m = get_model("liberty")
    base = CASES["liberty"][0]
    lo = m.simulate({**base, "thick": k - 1e-9})
    hi = m.simulate({**base, "thick": k + 1e-9})
    # Benford's fractional step goes from 1 to 2 layers as fraction 0 -> 1;
    # the 64/3 cubic correction makes the junction continuous to ~1e-3.
    assert np.abs(lo[0] - hi[0]).max() < 5e-3
    assert np.abs(lo[1] - hi[1]).max() < 5e-3


def test_liberty_no_real_root_region_is_nan():
    m = get_model("liberty")
    R, _ = m.simulate({**CASES["liberty"][0], "xu": 0.1, "baseline": 0.0001})
    assert np.isnan(R).any()


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


def test_liberty_absorber_order():
    assert ABSORBERS == ("c_factor", "w_factor", "element", "l_factor", "p_factor")
