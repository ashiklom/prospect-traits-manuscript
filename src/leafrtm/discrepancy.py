"""Learn shared RTM discrepancy from first-stage residuals.

Stage 1 fits a calibration pool of leaves with the ``ou`` error model. Its
posterior-mean residuals ``y - yhat`` (on each leaf's observed grid) are the
input here. For each (RTM, instrument class) we learn, on a 5-nm reference
grid over 400-2500 nm and jointly for the stacked [R; T] vector (so factors
carry R/T cross-correlation):

* ``mu``: mean bias (shared misfit structure)
* ``Lambda``: k factor loadings (probabilistic PCA with missing values)
* ``D``: per-band variance left after the factors
* ``ell``, ``nugget``: Matern-3/2 length scale and nugget of the
  leftover, from its empirical correlogram

The same correlogram fit, applied to stage-1 residuals standardized by their
heteroskedastic scale, gives the fixed correlation used by the ``ss`` model.

``map_to_obs`` interpolates a learned discrepancy onto a leaf's observed
bands and precomputes the Vecchia/Woodbury factors ``loglik_ld`` needs.
Cross-fitting: callers learn one discrepancy per held-out dataset fold, and
each leaf uses the one learned without its own fold.

Reference implementation: after commit 0cc1e51 (tag v0.3-pilot-all-errmodels) the
workflow fits hetero + OU errors only and does not run this module. It is kept,
with its tests, to document and reproduce the pilot's error-model comparison.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy.optimize import least_squares

from leafrtm import likelihood as L
from leafrtm.obsop import QCODE, ObsData

REF_GRID = np.arange(400.0, 2501.0, 5.0)
LD_QUANTITIES = ("R", "T")


@dataclass
class Discrepancy:
    quantities: tuple[str, ...]
    grid: np.ndarray
    mu: np.ndarray  # (nq * ngrid,)
    Lambda: np.ndarray  # (nq * ngrid, k)
    D: np.ndarray  # (nq * ngrid,)
    ell: float
    nugget: float
    n_leaves: int

    def save(self, path):
        d = asdict(self)
        d["quantities"] = np.array(self.quantities)
        np.savez(path, **d)

    @classmethod
    def load(cls, path):
        z = np.load(path)
        return cls(quantities=tuple(str(q) for q in z["quantities"]), grid=z["grid"],
                   mu=z["mu"], Lambda=z["Lambda"], D=z["D"], ell=float(z["ell"]),
                   nugget=float(z["nugget"]), n_leaves=int(z["n_leaves"]))


def to_ref(wavelength, values, grid=REF_GRID):
    """Interpolate a residual spectrum onto the reference grid (NaN outside
    the observed range, and across observed gaps wider than 20 nm)."""
    wavelength = np.asarray(wavelength, float)
    out = np.interp(grid, wavelength, values, left=np.nan, right=np.nan)
    j = np.clip(np.searchsorted(wavelength, grid), 1, len(wavelength) - 1)
    gap = wavelength[j] - wavelength[j - 1]
    out[(gap > 20) & ~np.isin(grid, wavelength)] = np.nan
    return out


def residual_matrix(leaves, quantities=LD_QUANTITIES, grid=REF_GRID):
    """Rows = leaves, columns = stacked quantities x grid (NaN where missing).

    ``leaves``: iterable of dicts {quantity: (wavelength, residual)}.
    """
    rows = []
    for leaf in leaves:
        parts = [to_ref(*leaf[q], grid) if q in leaf else np.full(len(grid), np.nan)
                 for q in quantities]
        rows.append(np.concatenate(parts))
    return np.array(rows)


def ppca_missing(X, k, n_iter=50):
    """Probabilistic PCA with missing values by iterative imputation.

    Returns mu, Lambda (p, k) and per-column residual variance D (p,)."""
    obs = np.isfinite(X)
    mu = np.nanmean(X, axis=0)
    mu = np.where(np.isfinite(mu), mu, 0.0)
    Xc = np.where(obs, X - mu, 0.0)
    for _ in range(n_iter):
        U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
        recon = (U[:, :k] * s[:k]) @ Vt[:k]
        Xc = np.where(obs, X - mu, recon)
    n = X.shape[0]
    U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
    ev = s**2 / max(n - 1, 1)
    noise = ev[k:].mean() if len(ev) > k else 0.0
    Lam = Vt[:k].T * np.sqrt(np.maximum(ev[:k] - noise, 0.0))
    resid = Xc - (U[:, :k] * s[:k]) @ Vt[:k]
    resid = np.where(obs, resid, np.nan)
    D = np.nanvar(resid, axis=0)
    floor = np.nanmedian(D) * 1e-2
    D = np.where(np.isfinite(D), np.maximum(D, floor), np.nanmax(D))
    return mu, Lam, D, resid


def fit_correlogram(series, max_lag_nm=150.0):
    """Fit Matern-3/2 (ell, nugget) to the pooled empirical correlogram of
    standardized residual series. ``series``: iterable of (wavelength, z) on
    regular grids (each normalized internally)."""
    lags, num, den = {}, {}, {}
    for wl, z in series:
        ok = np.isfinite(z)
        wl, z = np.asarray(wl)[ok], np.asarray(z)[ok]
        if len(z) < 10:
            continue
        step = np.median(np.diff(wl))
        z = (z - z.mean()) / (z.std() + 1e-300)
        for h in range(0, int(max_lag_nm / step) + 1):
            lag = round(h * step)
            num[lag] = num.get(lag, 0.0) + np.sum(z[: len(z) - h] * z[h:])
            den[lag] = den.get(lag, 0.0) + (len(z) - h)
            lags[lag] = lag
    lag = np.array(sorted(lags))
    rho = np.array([num[g] / den[g] for g in lag])

    def resid(p):
        ell, nug = np.exp(p[0]), 1 / (1 + np.exp(-p[1]))
        return L.matern32(lag[1:], ell, nug) - rho[1:]

    sol = least_squares(resid, x0=[np.log(30.0), -3.0])
    return float(np.exp(sol.x[0])), float(1 / (1 + np.exp(-sol.x[1]))), lag, rho


def learn(leaves, k=6, quantities=LD_QUANTITIES, grid=REF_GRID) -> Discrepancy:
    X = residual_matrix(leaves, quantities, grid)
    keep = np.isfinite(X).mean(axis=0) > 0.05
    if keep.sum() < k + 2 or X.shape[0] < k + 2:
        raise ValueError("Too few residuals to learn a discrepancy")
    mu, Lam, D, resid = ppca_missing(X, k)
    ng = len(grid)
    series = []
    for row in resid:
        for qi in range(len(quantities)):
            seg = row[qi * ng:(qi + 1) * ng] / np.sqrt(D[qi * ng:(qi + 1) * ng])
            series.append((grid, seg))
    ell, nugget, _, _ = fit_correlogram(series)
    return Discrepancy(tuple(quantities), grid, mu, Lam, D, ell, nugget, X.shape[0])


def map_to_obs(disc: Discrepancy, od: ObsData, m=30) -> dict:
    """Arrays for ``loglik_ld`` on this observation's bands."""
    n, ng = len(od), len(disc.grid)
    mu = np.zeros(n)
    D = np.ones(n)
    Lam = np.zeros((n, disc.Lambda.shape[1]))
    for qi, q in enumerate(disc.quantities):
        sel = od.qcode == QCODE[q]
        if not sel.any():
            continue
        sl = slice(qi * ng, (qi + 1) * ng)
        x = od.wavelength[sel]
        mu[sel] = np.interp(x, disc.grid, disc.mu[sl])
        D[sel] = np.interp(x, disc.grid, disc.D[sl])
        for j in range(Lam.shape[1]):
            Lam[sel, j] = np.interp(x, disc.grid, disc.Lambda[sl, j])
    other = ~np.isin(od.qcode, [QCODE[q] for q in disc.quantities])
    if other.any():
        raise ValueError("Discrepancy does not cover all observed quantities")
    sqrtD = np.sqrt(D)
    nb, b, S = vecchia_blocks(od, disc.ell, disc.nugget, m=m)
    W, Minv, logdetM = L.ld_factors(Lam / sqrtD[:, None], nb, b, S)
    return {"mu": mu, "sqrtD": sqrtD, "nb": nb, "b": b, "S": S, "W": W, "Minv": Minv,
            "logdetM": np.float64(logdetM)}


def vecchia_blocks(od: ObsData, ell, nugget, m=30):
    """Vecchia factors for every block of an observation (flat indexing)."""
    corr = lambda d: L.matern32(d, ell, nugget)  # noqa: E731
    nb = np.zeros((len(od), m), dtype=np.int64)
    b = np.zeros((len(od), m))
    S = np.ones(len(od))
    for bi in np.unique(od.block):
        idx = np.flatnonzero(od.block == bi)
        n_, b_, S_ = L.vecchia(od.wavelength[idx], od.start[idx], corr, m=m)
        nb[idx] = idx[n_]
        b[idx] = b_
        S[idx] = S_
    return nb, b, S
