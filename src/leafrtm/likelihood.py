"""Error models (likelihoods) for spectral residuals.

Residuals ``e = y - yhat`` are modeled per block (see ``leafrtm.obsop``). All
likelihoods are exact or near-exact and cost O(n) (O(n k) for the low-rank
discrepancy), and all work on irregular, gappy wavelength grids.

``ou``   Heteroskedastic Ornstein-Uhlenbeck (continuous AR(1)) errors:
         ``e_i = s_i z_i``, ``s_i = sigma0 + sigma1 |yhat_i|``, and z a unit-variance
         OU process with correlation ``exp(-|dlam| / ell)``. sigma0, sigma1 and ell
         are free per quantity. This generalizes the heteroskedastic AR(1)
         model of prospect-inversion-comparison to irregular grids.

``ss``   Heteroskedastic smooth state-space errors: as ``ou`` but z is a
         Matern-3/2 process (plus a small nugget) whose length scale and nugget
         are fixed (learned from first-stage residuals). With the correlation
         fixed, the exact Gaussian density is evaluated through a precomputed
         Vecchia factorization: each band conditions on its ``m`` predecessors
         in the same segment. With m = 30 the log density is within ~0.01 nats
         of the dense density for residuals consistent with the model (tested);
         the nugget weakens Vecchia's screening, so m must not be too small.

``ld``   Learned shared discrepancy: ``e - mu_g ~ N(0, kappa^2 (D^1/2 C D^1/2 +
         Lambda Lambda^T))``. The mean bias mu_g, per-band scale D_g, factor
         loadings Lambda (k factors), and Matern correlation C are learned
         across leaves for each (RTM, quantity, instrument class)
         (``leafrtm.discrepancy``); only the per-quantity scale kappa is free per
         leaf. It is evaluated with Vecchia for B = D^1/2 C D^1/2 and Woodbury
         for the low-rank part.

Correlated R/T errors (``rt``):
  * for ``ou``/``ss``, a separable bivariate process: T's standardized
    residual is ``rho z_R + sqrt(1 - rho^2) xi`` on the same grid, with xi
    independent of z_R and with R's correlation structure;
  * for ``ld``, factors learned on the stacked [R; T] residual vector carry
    the cross-correlation.
"""

from __future__ import annotations

import numpy as np
import pytensor.tensor as pt

LOG2PI = float(np.log(2 * np.pi))


# ---------------------------------------------------------------------------
# Correlation functions and Vecchia precomputation (NumPy)
# ---------------------------------------------------------------------------


def matern32(d, ell, nugget=0.0):
    """Matern-3/2 correlation with a nugget (variance fraction ``nugget``)."""
    d = np.abs(np.asarray(d, dtype=float))
    a = np.sqrt(3.0) * d / ell
    c = (1 - nugget) * (1 + a) * np.exp(-a)
    return np.where(d == 0, 1.0, c)


def segment_ids(start):
    return np.cumsum(np.asarray(start, dtype=bool)) - 1


def dense_corr(wavelength, start, corr):
    """Dense correlation matrix (zero across segments); for tests and learning."""
    seg = segment_ids(start)
    C = corr(wavelength[:, None] - wavelength[None, :])
    C[seg[:, None] != seg[None, :]] = 0.0
    return C


def vecchia(wavelength, start, corr, m=30):
    """Vecchia factors for a fixed correlation function.

    Returns (nb, b, S): for band i, its predecessors nb[i] (padded with i,
    coefficient 0), regression coefficients b[i] and conditional variance S[i],
    so that ``(z_i - b_i . z[nb_i]) / sqrt(S_i)`` are iid N(0, 1).
    Bands in different segments are independent.
    """
    wl = np.asarray(wavelength, dtype=float)
    seg = segment_ids(start)
    n = len(wl)
    nb = np.tile(np.arange(n)[:, None], (1, m))
    b = np.zeros((n, m))
    S = np.ones(n)
    seg_first = np.r_[0, np.flatnonzero(np.diff(seg)) + 1]
    first = seg_first[seg]
    for i in range(n):
        J = np.arange(max(i - m, first[i]), i)
        if len(J) == 0:
            continue
        Cjj = corr(wl[J, None] - wl[None, J])
        cji = corr(wl[J] - wl[i])
        bi = np.linalg.solve(Cjj, cji)
        nb[i, : len(J)] = J
        b[i, : len(J)] = bi
        S[i] = 1.0 - cji @ bi
    return nb, b, S


# ---------------------------------------------------------------------------
# PyTensor building blocks
# ---------------------------------------------------------------------------


def hetero_scale(yhat, qcode, sigma0, sigma1):
    """s_i = sigma0[q_i] + sigma1[q_i] |yhat_i| (one sigma pair per quantity)."""
    return sigma0[qcode] + sigma1[qcode] * pt.abs(yhat)


def rt_adjust(z, pair, rho):
    """Separable correlated R/T: replace paired T values by their innovation
    given R, (z_T - rho z_R) / sqrt(1 - rho^2). Returns (z_adj, logdet term)."""
    paired = pair >= 0
    src = pt.switch(paired, pair, pt.arange(z.shape[0]))
    scale = pt.sqrt(1 - rho**2)
    z_adj = pt.switch(paired, (z - rho * z[src]) / scale, z)
    return z_adj, pt.sum(paired) * pt.log(scale)


def ou_loglik(z, dlam, start, ell_i):
    """Exact log density of a unit-variance OU process sampled at irregular
    spacings ``dlam`` (restarting at ``start``); ``ell_i`` is the length scale
    for each band."""
    x = -dlam / ell_i
    phi = pt.switch(start, 0.0, pt.exp(x))
    one_m_phi2 = pt.switch(start, 1.0, -pt.expm1(2 * x))
    log_one_m_phi2 = pt.switch(start, 0.0, pt.log1mexp(2 * x))
    z_prev = pt.concatenate([pt.zeros((1,)), z[:-1]])
    innov = z - phi * z_prev
    return -0.5 * pt.sum(LOG2PI + log_one_m_phi2 + innov**2 / one_m_phi2)


def vecchia_whiten(z, nb, b, S):
    return (z - pt.sum(b * z[nb], axis=1)) / pt.sqrt(S)


def vecchia_loglik(z, nb, b, S):
    u = vecchia_whiten(z, nb, b, S)
    return -0.5 * pt.sum(LOG2PI + pt.log(S) + u**2)


# ---------------------------------------------------------------------------
# Full error models
# ---------------------------------------------------------------------------


def loglik_ou(e, yhat, d, sigma0, sigma1, ell, rho=None):
    """d: dict with qcode, dlam, start, pair. sigma0/sigma1/ell indexed by qcode.
    With ``rho``, T uses R's length scale (separable bivariate process)."""
    s = hetero_scale(yhat, d["qcode"], sigma0, sigma1)
    z = e / s
    logdet = pt.sum(pt.log(s))
    ell_i = ell[d["qcode"]]
    if rho is not None:
        z, ld_rt = rt_adjust(z, d["pair"], rho)
        logdet += ld_rt
        ell_i = pt.switch(d["pair"] >= 0, ell[0], ell_i)
    return ou_loglik(z, d["dlam"], d["start"], ell_i) - logdet


def loglik_ss(e, yhat, d, sigma0, sigma1, rho=None):
    """d additionally holds the fixed Vecchia factors nb, b, S."""
    s = hetero_scale(yhat, d["qcode"], sigma0, sigma1)
    z = e / s
    logdet = pt.sum(pt.log(s))
    if rho is not None:
        z, ld_rt = rt_adjust(z, d["pair"], rho)
        logdet += ld_rt
    return vecchia_loglik(z, d["nb"], d["b"], d["S"]) - logdet


def loglik_ld(e, d, kappa):
    """d holds the learned discrepancy on the observed grid: mu, sqrtD, the
    Vecchia factors (nb, b, S) of C, W = B^-1/2 Lambda (n, k), Minv =
    (I + W^T W)^-1 and logdetM. kappa is indexed by qcode."""
    kq = kappa[d["qcode"]]
    r = (e - d["mu"]) / (kq * d["sqrtD"])
    u = vecchia_whiten(r, d["nb"], d["b"], d["S"])
    c = pt.dot(d["W"].T, u)
    quad = pt.sum(u**2) - pt.dot(c, pt.dot(d["Minv"], c))
    logdet = (2 * pt.sum(pt.log(kq * d["sqrtD"])) + pt.sum(pt.log(d["S"])) + d["logdetM"])
    return -0.5 * (e.shape[0] * LOG2PI + logdet + quad)


def ld_factors(W_raw, nb, b, S):
    """Precompute W = C^-1/2 (Lambda / sqrtD) and (I + W^T W)^-1 for ``loglik_ld``.

    ``W_raw`` is Lambda / sqrtD on the observed grid (columns = factors);
    whitening applies the Vecchia transform column by column."""
    W = (W_raw - np.einsum("nm,nmk->nk", b, W_raw[nb])) / np.sqrt(S)[:, None]
    M = np.eye(W.shape[1]) + W.T @ W
    sign, logdetM = np.linalg.slogdet(M)
    assert sign > 0
    return W, np.linalg.inv(M), logdetM
