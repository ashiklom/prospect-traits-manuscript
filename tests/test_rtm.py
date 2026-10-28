"""Forward-model tests for all leaf RTMs."""

import numpy as np
import pytensor
import pytensor.tensor as pt
import pytest

from leafrtm.rtm import MODELS, get_model
from leafrtm.rtm.e1 import e1_fast_np
from leafrtm.rtm.gpm import tav
from leafrtm.rtm.liberty import ABSORBERS, _cell, _slab_m_np
from liberty_rrtm import liberty_rrtm
from scipy.special import exp1
from sip_matlab import sip_model

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
    # the second case has no real R_inf in LIBERTY as published
    "liberty_fixed": [dict(xu=0.045, thick=1.6, baseline=0.0006, element=2, c_factor=200, w_factor=100,
                           l_factor=40, p_factor=1),
                      dict(xu=0.1, thick=2.5, baseline=0.0003, element=1, c_factor=200, w_factor=100,
                           l_factor=40, p_factor=1)],
    "leafsip": [dict(Cab=30, Car=10, Canth=1, Cbrown=0, Cw=0.015, Cm=0.009),
                dict(Cab=60, Car=15, Canth=3, Cbrown=0.3, Cw=0.03, Cm=0.02)],
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


@pytest.mark.parametrize("case", CASES["leafsip"])
def test_leafsip_matches_matlab_transliteration(case):
    m = get_model("leafsip")
    R, T = m.simulate(case)
    Rm, Tm = sip_model(m.K, *(case[c] for c in ("Cab", "Car", "Canth", "Cbrown", "Cw", "Cm")))
    np.testing.assert_allclose(R, Rm, atol=1e-12)
    np.testing.assert_allclose(T, Tm, atol=1e-12)


def _interface(m):
    """The oracle with this port's (corrected) interface reflectances."""
    return dict(me_override=m.me, mi_override=m.mi)


def _fresnel(theta, n_from, n_to):
    """Fresnel reflectance of unpolarized light at incidence angle theta (in
    the medium of index n_from); 1 beyond the critical angle."""
    s = n_from / n_to * np.sin(theta)
    out = np.arcsin(np.minimum(s, 1))
    p, d = theta + out, theta - out
    r = 0.5 * (np.sin(d) ** 2 / np.sin(p) ** 2 + np.tan(d) ** 2 / np.tan(p) ** 2)
    return np.where(s >= 1, 1.0, r)


@pytest.mark.parametrize("i", [0, 80, 320, 419])
def test_liberty_interface_reflectances_are_dawson_eqs_1_2(i):
    # Eq. 1 (with the factor 2 that the printed equation omits) and Eq. 2, by
    # Gauss-Legendre quadrature, with the Fresnel reflectance of light from the
    # air (me) and from inside the cell (mi)
    x, w = np.polynomial.legendre.leggauss(400)

    def integral(f, lo, hi):
        t = (hi - lo) / 2 * x + (hi + lo) / 2
        return (hi - lo) / 2 * np.sum(w * f(t))

    n = 1.4891 - 0.0005 * i
    crit = np.arcsin(1 / n)
    me = integral(lambda t: 2 * _fresnel(t, 1, n) * np.sin(t) * np.cos(t), 0, np.pi / 2)
    mi = (1 - np.sin(crit) ** 2) + integral(lambda t: 2 * _fresnel(t, n, 1) * np.sin(t) * np.cos(t),
                                            0, crit)
    np.testing.assert_allclose([me, mi], [1 - tav(90, n), 1 - tav(90, n) / n**2], atol=1e-10)
    m = get_model("liberty")
    np.testing.assert_allclose([m.me[i], m.mi[i]], [me, mi], atol=1e-10)


@pytest.mark.parametrize("case", CASES["liberty"])
def test_liberty_matches_converged_rrtm(case):
    m = get_model("liberty")
    R, T, Rinf = m.simulate(case, return_rinf=True)
    Ro, To, Rinfo = liberty_rrtm(m.K, m.fixed["D"], **case, n_iter=20000, **_interface(m))
    np.testing.assert_allclose(R, Ro, atol=1e-11)
    np.testing.assert_allclose(T, To, atol=1e-11)
    np.testing.assert_allclose(Rinf, Rinfo, atol=1e-11)
    # rrtm's 50 iterations are close to, but not exactly at, the fixed point
    Ro50, To50, _ = liberty_rrtm(m.K, m.fixed["D"], **case, **_interface(m))
    assert np.abs(R - Ro50).max() < 1e-4


@pytest.mark.parametrize("thick", [1.2, 2.0, 2.7, 3.0, 4.0, 6.5])
def test_liberty_integer_layers_match_loop(thick):
    m = get_model("liberty")
    case = {**CASES["liberty"][0], "thick": thick}
    R, T = m.simulate(case)
    Ro, To, _ = liberty_rrtm(m.K, m.fixed["D"], **case, n_iter=20000, **_interface(m))
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


def test_liberty_fixed_rinf_is_the_sum_of_dawson_eq_15():
    # R_inf must reproduce itself as the R under the surface layer in Eq. 15,
    # summed term by term; Eq. 16 (LIBERTY as published) does not
    m = get_model("liberty_fixed")
    rng = np.random.default_rng(0)
    for _ in range(50):
        i = rng.integers(420)
        me, mi = m.me[i], m.mi[i]
        x, tau = _cell(mi, rng.uniform(0.005, 0.5), _slab_m_np(np.exp(rng.uniform(-8, 2))))
        rd, tb, R = m._layer(me, x, tau, np)
        q = (1 - me) * tau * R / (1 - me * R)
        eq15 = 2 * x * me + sum(x * (1 - x) ** j * (1 - 2 * x * me) * tau * q**j for j in range(5000))
        assert eq15 == pytest.approx(R, rel=1e-12)
        assert tb**2 == pytest.approx((R - rd) * (1 - R * rd) / R, rel=1e-9, abs=1e-15)  # Eq. 19
        with np.errstate(invalid="ignore"):
            R16 = get_model("liberty")._layer(me, x, tau, np)[2]
        if np.isfinite(R16):
            assert not np.isclose(R16, R, rtol=1e-6, atol=0)


def test_liberty_fixed_defined_everywhere_and_conserves_energy():
    m, pub = get_model("liberty_fixed"), get_model("liberty")
    zero = dict(element=0, c_factor=0, w_factor=0, l_factor=0, p_factor=0)
    n_pub_undefined = 0
    for xu in (0.005, 0.03, 0.06, 0.1, 0.3, 0.5):
        for baseline in (1e-8, 1e-6, 1e-4, 1e-3, 1e-2):
            for thick in (1.0, 1.6, 3.4, 10.0):
                case = {**CASES["liberty"][0], "xu": xu, "baseline": baseline, "thick": thick}
                R, T = m.simulate(case)
                assert np.all(np.isfinite(R) & np.isfinite(T)), case
                assert np.all((R > 0) & (T > 0) & (R + T <= 1 + 1e-12)), case
                n_pub_undefined += np.isnan(pub.simulate(case)[0]).any()
        # no absorption: no light is lost (R and T themselves are 0/0 there)
        _, _, Rinf = m.simulate({**zero, "xu": xu, "baseline": 0.0, "thick": 2.3}, return_rinf=True)
        np.testing.assert_allclose(Rinf, 1, atol=1e-12)
        R, T = m.simulate({**zero, "xu": xu, "baseline": 1e-8, "thick": 2.3})
        np.testing.assert_allclose(R + T, 1, atol=1e-5)
    assert n_pub_undefined > 0
    assert m.validity_margin(CASES["liberty"][0]) is None


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
