import numpy as np
import pytensor
import pytensor.tensor as pt

from leafrtm import likelihood as L
from leafrtm.discrepancy import REF_GRID, Discrepancy, learn, map_to_obs
from leafrtm.obsop import Block, make_obsdata


def _synthetic(n_leaves=300, seed=0):
    rng = np.random.default_rng(seed)
    wl = np.arange(400.0, 2501.0, 2.0)
    g = (wl - 400) / 2100
    mu = {"R": 0.01 * np.sin(6 * g), "T": -0.008 * np.cos(4 * g)}
    f1 = np.exp(-((wl - 680) / 40) ** 2)  # red-edge-like factor shared by R and T
    f2 = np.exp(-((wl - 1450) / 60) ** 2)
    lam = {"R": np.c_[0.01 * f1, 0.006 * f2], "T": np.c_[-0.008 * f1, 0.004 * f2]}
    C = L.matern32(wl[:, None] - wl[None, :], 40.0, 0.05)
    Lc = np.linalg.cholesky(C + 1e-10 * np.eye(len(wl)))
    leaves = []
    for _ in range(n_leaves):
        z = rng.normal(size=2)
        leaf = {}
        for q in ("R", "T"):
            e = mu[q] + lam[q] @ z + 0.002 * (Lc @ rng.normal(size=len(wl)))
            leaf[q] = (wl, e)
        if rng.uniform() < 0.3:
            del leaf["T"]  # some leaves are R-only
        leaves.append(leaf)
    return leaves, mu, lam, wl


def test_learn_recovers_structure(tmp_path):
    leaves, mu, lam, wl = _synthetic()
    disc = learn(leaves, k=2)
    ng = len(REF_GRID)
    i = ((REF_GRID - 400) / 2).astype(int)
    assert np.corrcoef(disc.mu[:ng], mu["R"][i])[0, 1] > 0.95
    assert np.corrcoef(disc.mu[ng:], mu["T"][i])[0, 1] > 0.95
    true = np.r_[lam["R"][i], lam["T"][i]]
    # learned factor subspace spans the true loadings
    Q, _ = np.linalg.qr(disc.Lambda)
    resid = true - Q @ (Q.T @ true)
    assert np.linalg.norm(resid) / np.linalg.norm(true) < 0.15
    assert 20 < disc.ell < 80
    disc.save(tmp_path / "d.npz")
    d2 = Discrepancy.load(tmp_path / "d.npz")
    np.testing.assert_array_equal(d2.Lambda, disc.Lambda)


def test_map_to_obs_gives_finite_loglik():
    leaves, *_ = _synthetic(80, seed=1)
    disc = learn(leaves, k=2)
    wl = np.arange(400.0, 2501.0, 1.0)
    od = make_obsdata([Block("R", wl, np.full(len(wl), 0.3)), Block("T", wl, np.full(len(wl), 0.3))], wl)
    arrs = map_to_obs(disc, od)
    d = {k: pt.as_tensor_variable(v) for k, v in arrs.items()}
    d["qcode"] = pt.as_tensor_variable(od.qcode)
    E, K = pt.vector(), pt.vector()
    f = pytensor.function([E, K], L.loglik_ld(E, d, K))
    val = f(arrs["mu"], np.array([1.0, 1.0, 1, 1]))
    assert np.isfinite(val)
