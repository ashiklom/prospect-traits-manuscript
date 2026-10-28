"""Generalized plate model (Allen et al. 1969; Stokes 1862) and interface
transmissivity (Stern 1964; Allen 1973).

Vendored from ashiklom/prospect-inversion-comparison @ 84504ed
(``python/prospect_np.py``, ``python/prospect_pt.py``), which in turn ports
``gpm.R``/``tav_abs`` from ashiklom/rrtm (MIT).
"""

from __future__ import annotations

import numpy as np
import pytensor.tensor as pt

from leafrtm.rtm.e1 import E1_NP, E1_PT


def tav(theta_deg, n):
    """Average transmissivity of a dielectric plane surface.

    ``theta_deg`` is the maximum incidence angle (degrees), ``n`` the
    refractive index (array).
    """
    n = np.asarray(n, dtype=float)
    theta = np.deg2rad(theta_deg)
    if theta_deg == 0:
        return 4 * n / (n + 1) ** 2
    r2 = n**2
    rp = r2 + 1
    rm = r2 - 1
    a = (n + 1) ** 2 / 2
    k = -((r2 - 1) ** 2) / 4
    sa = np.sin(theta)
    b1 = 0.0 if theta_deg == 90 else np.sqrt((sa**2 - rp / 2) ** 2 + k)
    b2 = sa**2 - rp / 2
    b = b1 - b2
    ts = (k**2 / (6 * b**3) + k / b - b / 2) - (k**2 / (6 * a**3) + k / a - a / 2)
    tp1 = -2 * r2 * (b - a) / rp**2
    tp2 = -2 * r2 * rp * np.log(b / a) / rm**2
    tp3 = r2 * (1 / b - 1 / a) / 2
    tp4 = (16 * r2**2 * (r2**2 + 1) * np.log((2 * rp * b - rm**2) / (2 * rp * a - rm**2))
           / (rp**3 * rm**2))
    tp5 = 16 * r2**3 * (1 / (2 * rp * b - rm**2) - 1 / (2 * rp * a - rm**2)) / rp**3
    tp = tp1 + tp2 + tp3 + tp4 + tp5
    return (ts + tp) / (2 * sa**2)


def interface_terms(n):
    """``talf = tav(40, n)``, ``t12 = tav(90, n)``, ``t21 = t12 / n**2``."""
    t12 = tav(90, n)
    return tav(40, n), t12, t12 / np.asarray(n) ** 2


# ---------------------------------------------------------------------------
# NumPy reference
# ---------------------------------------------------------------------------


def transmissivity_np(k, e1="fast"):
    """Transmissivity of a compact layer with absorption coefficient k."""
    k = np.asarray(k, dtype=float)
    kpos = np.where(k > 0, k, 1.0)
    tau = (1 - kpos) * np.exp(-kpos) + kpos**2 * E1_NP[e1](kpos)
    tau = np.where(k > 0, tau, 1.0)
    return np.maximum(tau, 0.0)


def stokes_np(r, t, n):
    """Reflectance and transmittance of ``n`` identical layers (r, t)."""
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        D = np.sqrt(np.maximum((1 + r + t) * (1 + r - t) * (1 - r + t) * (1 - r - t), 0))
        va = (1 + r**2 - t**2 + D) / (2 * r)
        vb = (1 - r**2 + t**2 + D) / (2 * t)
        vbNN = vb**n
        denomx = va**2 * vbNN**2 - 1
        Rsub = va * (vbNN**2 - 1) / denomx
        Tsub = vbNN * (va**2 - 1) / denomx
    # No absorption (r + t = 1): Kubelka-Munk-like limit
    noabs = (r + t) >= 1
    Tsub_na = t / (t + (1 - t) * n)
    Rsub = np.where(noabs, 1 - Tsub_na, Rsub)
    Tsub = np.where(noabs, Tsub_na, Tsub)
    # Total absorption
    opaque = (t <= 0) | (r <= 0)
    Rsub = np.where(opaque & ~noabs, 0.0, Rsub)
    Tsub = np.where(opaque & ~noabs, 0.0, Tsub)
    return Rsub, Tsub


def gpm_np(k, N, talf, t12, t21, e1="fast"):
    """Generalized plate model: (R, T) of N plates with absorption k."""
    tau = transmissivity_np(k, e1)
    ralf = 1 - talf
    r12 = 1 - t12
    r21 = 1 - t21
    denom = 1 - r21**2 * tau**2
    Ta = talf * tau * t21 / denom
    Ra = ralf + r21 * tau * Ta
    t = t12 * tau * t21 / denom
    r = r12 + r21 * tau * t
    Rsub, Tsub = stokes_np(r, t, N - 1)
    denomy = 1 - Rsub * r
    return Ra + Ta * Rsub * t / denomy, Ta * Tsub / denomy


# ---------------------------------------------------------------------------
# PyTensor
# ---------------------------------------------------------------------------


def gpm_pt(k, N, talf, t12, t21, e1="fast"):
    """Generalized plate model in PyTensor; returns (R, T). Assumes k > 0
    (true whenever Cw or Cm > 0)."""
    tau = pt.maximum((1 - k) * pt.exp(-k) + k**2 * E1_PT[e1](k), 0.0)
    ralf = 1 - talf
    r12 = 1 - t12
    r21 = 1 - t21
    denom = 1 - r21**2 * tau**2
    Ta = talf * tau * t21 / denom
    Ra = ralf + r21 * tau * Ta
    t = t12 * tau * t21 / denom
    r = r12 + r21 * tau * t
    D = pt.sqrt((1 + r + t) * (1 + r - t) * (1 - r + t) * (1 - r - t))
    va = (1 + r**2 - t**2 + D) / (2 * r)
    vb = (1 - r**2 + t**2 + D) / (2 * t)
    vbNN = vb ** (N - 1)
    vbNN2 = vbNN**2
    denomx = va**2 * vbNN2 - 1
    Rsub = va * (vbNN2 - 1) / denomx
    Tsub = vbNN * (va**2 - 1) / denomx
    denomy = 1 - Rsub * r
    return Ra + Ta * Rsub * t / denomy, Ta * Tsub / denomy
