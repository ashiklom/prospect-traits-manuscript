"""Likelihoods vs. dense multivariate normal densities on irregular grids."""

import numpy as np
import pytensor
import pytensor.tensor as pt
import pytest
from scipy.stats import multivariate_normal

from leafrtm import likelihood as L
from leafrtm.obsop import Block, make_obsdata

rng = np.random.default_rng(42)


def _grid():
    # irregular, gappy grid with a detector segment boundary
    wl = np.r_[np.arange(400, 700, 1.0), np.arange(705, 1000, 2.0), np.arange(1000, 1300, 1.5),
               np.arange(1500, 1800, 3.0)]
    start = np.zeros(len(wl), bool)
    start[0] = True
    start[np.searchsorted(wl, 1000)] = True  # detector join
    return wl, start


def _ou_corr(wl, start, ell):
    return L.dense_corr(wl, start, lambda d: np.exp(-np.abs(d) / ell))


def _obs(two_blocks=True):
    wl, start = _grid()
    model_wl = np.arange(400, 2501, 1.0)
    blocks = [Block("R", wl, rng.uniform(0.05, 0.5, len(wl)), segment_start=start)]
    if two_blocks:
        blocks.append(Block("T", wl, rng.uniform(0.05, 0.5, len(wl)), segment_start=start))
    return make_obsdata(blocks, model_wl)


def _fn(expr, inputs):
    return pytensor.function(inputs, expr)


@pytest.mark.parametrize("with_rho", [False, True])
def test_ou_matches_dense(with_rho):
    od = _obs()
    yhat = rng.uniform(0.05, 0.5, len(od))
    e = rng.normal(0, 0.01, len(od))
    s0, s1, ell, rho = np.array([0.004, 0.006, 0, 0]), np.array([0.05, 0.02, 0, 0]), np.array([35.0, 80.0, 1, 1]), 0.6
    d = {k: pt.as_tensor_variable(getattr(od, k)) for k in ("qcode", "dlam", "start", "pair")}
    E, Y, S0, S1, ELL, RHO = pt.vector(), pt.vector(), pt.vector(), pt.vector(), pt.vector(), pt.scalar()
    expr = L.loglik_ou(E, Y, d, S0, S1, ELL, rho=RHO if with_rho else None)
    f = pytensor.function([E, Y, S0, S1, ELL, RHO], expr, on_unused_input="ignore")
    got = f(e, yhat, s0, s1, ell, rho)

    s = s0[od.qcode] + s1[od.qcode] * yhat
    isR, isT = od.qcode == 0, od.qcode == 1
    C = np.zeros((len(od), len(od)))
    if with_rho:
        CR = _ou_corr(od.wavelength[isR], od.start[isR], ell[0])
        C[np.ix_(isR, isR)] = CR
        C[np.ix_(isT, isT)] = CR
        C[np.ix_(isR, isT)] = rho * CR
        C[np.ix_(isT, isR)] = rho * CR
    else:
        C[np.ix_(isR, isR)] = _ou_corr(od.wavelength[isR], od.start[isR], ell[0])
        C[np.ix_(isT, isT)] = _ou_corr(od.wavelength[isT], od.start[isT], ell[1])
    want = multivariate_normal(np.zeros(len(od)), s[:, None] * C * s[None, :]).logpdf(e)
    assert got == pytest.approx(want, rel=1e-9)


def _vecchia_for(od, ell, nugget, m=30):
    corr = lambda dd: L.matern32(dd, ell, nugget)  # noqa: E731
    nb, b, S = [], [], []
    off = 0
    for bi in np.unique(od.block):
        sel = od.block == bi
        n_, b_, S_ = L.vecchia(od.wavelength[sel], od.start[sel], corr, m=m)
        nb.append(n_ + off)
        b.append(b_)
        S.append(S_)
        off += sel.sum()
    return np.concatenate(nb), np.concatenate(b), np.concatenate(S), corr


@pytest.mark.parametrize("with_rho", [False, True])
def test_ss_matches_dense(with_rho):
    od = _obs()
    nb, b, S, corr = _vecchia_for(od, ell=60.0, nugget=0.01)
    yhat = rng.uniform(0.05, 0.5, len(od))
    s0, s1, rho = np.array([0.004, 0.006, 0, 0]), np.array([0.05, 0.02, 0, 0]), -0.4
    d = {"qcode": od.qcode, "pair": od.pair, "nb": nb, "b": b, "S": S}
    d = {k: pt.as_tensor_variable(v) for k, v in d.items()}
    E, Y, S0, S1, RHO = pt.vector(), pt.vector(), pt.vector(), pt.vector(), pt.scalar()
    f = pytensor.function([E, Y, S0, S1, RHO],
                          L.loglik_ss(E, Y, d, S0, S1, rho=RHO if with_rho else None),
                          on_unused_input="ignore")
    s = s0[od.qcode] + s1[od.qcode] * yhat
    isR, isT = od.qcode == 0, od.qcode == 1
    CR = L.dense_corr(od.wavelength[isR], od.start[isR], corr)
    C = np.zeros((len(od), len(od)))
    C[np.ix_(isR, isR)] = CR
    C[np.ix_(isT, isT)] = CR
    if with_rho:
        C[np.ix_(isR, isT)] = rho * CR
        C[np.ix_(isT, isR)] = rho * CR
    Sigma = s[:, None] * C * s[None, :]
    e = np.linalg.cholesky(Sigma) @ rng.normal(size=len(od))
    got = f(e, yhat, s0, s1, rho)
    want = multivariate_normal(np.zeros(len(od)), Sigma).logpdf(e)
    assert got == pytest.approx(want, abs=0.1)


def test_ld_matches_dense():
    od = _obs()
    nb, b, S, corr = _vecchia_for(od, ell=40.0, nugget=0.02)
    n, k = len(od), 4
    mu = rng.normal(0, 0.003, n)
    sqrtD = rng.uniform(0.002, 0.01, n)
    Lam = rng.normal(0, 0.01, (n, k))  # spans R and T jointly (correlated R/T)
    W, Minv, logdetM = L.ld_factors(Lam / sqrtD[:, None], nb, b, S)
    kappa = np.array([1.3, 0.8, 1, 1])
    d = {"qcode": od.qcode, "mu": mu, "sqrtD": sqrtD, "nb": nb, "b": b, "S": S, "W": W,
         "Minv": Minv, "logdetM": logdetM}
    d = {k_: pt.as_tensor_variable(v) for k_, v in d.items()}
    E, K = pt.vector(), pt.vector()
    f = pytensor.function([E, K], L.loglik_ld(E, d, K))
    C = np.zeros((n, n))
    for q in (0, 1):
        sel = od.qcode == q
        C[np.ix_(sel, sel)] = L.dense_corr(od.wavelength[sel], od.start[sel], corr)
    kq = kappa[od.qcode]
    Sigma = kq[:, None] * (sqrtD[:, None] * C * sqrtD[None, :] + Lam @ Lam.T) * kq[None, :]
    e = mu + np.linalg.cholesky(Sigma) @ rng.normal(size=n)
    want = multivariate_normal(mu, Sigma).logpdf(e)
    assert f(e, kappa) == pytest.approx(want, abs=0.1)


def test_vecchia_is_exact_for_ou():
    # OU is Markov: one neighbor reproduces the exact density
    wl, start = _grid()
    ell = 50.0
    corr = lambda d: np.exp(-np.abs(d) / ell)  # noqa: E731
    nb, b, S = L.vecchia(wl, start, corr, m=1)
    z = rng.normal(size=len(wl))
    u = (z - np.sum(b * z[nb], 1)) / np.sqrt(S)
    got = -0.5 * np.sum(L.LOG2PI + np.log(S) + u**2)
    want = multivariate_normal(np.zeros(len(wl)), L.dense_corr(wl, start, corr)).logpdf(z)
    assert got == pytest.approx(want, rel=1e-10)


def test_partial_pairing_is_dropped():
    wl, start = _grid()
    model_wl = np.arange(400, 2501, 1.0)
    od = make_obsdata([Block("R", wl, np.full(len(wl), 0.3), segment_start=start),
                       Block("T", wl[::2], np.full(len(wl[::2]), 0.3))], model_wl)
    assert np.all(od.pair[od.qcode == 1] >= 0)  # T grid is a subset of R's grid
    od = make_obsdata([Block("R", wl[::2], np.full(len(wl[::2]), 0.3)),
                       Block("T", wl, np.full(len(wl), 0.3), segment_start=start)], model_wl)
    assert np.all(od.pair == -1)
