"""LIBERTY (Leaf Incorporating Biochemistry Exhibiting Reflectance and
Transmittance Yields; Dawson, Curran & Plummer 1998, RSE 65:50-60,
doi:10.1016/S0034-4257(98)00007-8).

Port of ``liberty()`` in ashiklom/rrtm (``R/liberty.R``, MIT), with that
package's coefficients (``data-raw/liberty_raw.dat``, Dawson et al. 1998
values; 400-2495 nm at 5 nm -- rrtm uses the first 420 of the 421 rows, and so
do we). Behavior-preserving refactors relative to rrtm, each checked against a
literal transliteration of the R loop (``tests/reference/liberty_rrtm.py``):

* The external/internal interface reflectances ``me`` and ``mi`` depend only on
  wavelength (through the refractive index), so they are precomputed.
* The infinite reflectance R_inf solves ``a R^2 + b R + c = 0``; rrtm finds it
  with 50 fixed-point iterations from 0.5. We take the attracting (smaller)
  root in closed, cancellation-free form, ``R = -2c / (b + sqrt(b^2 - 4ac))``.
  This is the *converged* fixed point; rrtm's 50 iterations are not fully
  converged where |dg/dR| is close to 1 (up to ~3e-6 at default ``xu``).
* The loop that stacks ``floor(thick) - 1`` identical layers is replaced by the
  closed-form Stokes N-layer solution (identical for integer layer counts).
* The 0/0 limit of the slab transmission ``M(coeff)`` as coeff -> 0 is handled
  with a series expansion.

Notes for inference: ``D`` (cell diameter) multiplies every absorption term,
so it is exactly non-identifiable from the concentrations; it is held fixed
(``fixed``) at Dawson's typical value. ``thick`` must be >= 1 (the model jumps
from 2 to 1 layers at thick = 1). For large ``xu`` with weak absorption the
quadratic has no real root (rrtm diverges); the PyTensor version returns NaN
there, which the sampler rejects.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from importlib.resources import files

import numpy as np
import pandas as pd
import pytensor.tensor as pt

from leafrtm.rtm.base import RTM

NWL = 420
ABSORBERS = ("c_factor", "w_factor", "element", "l_factor", "p_factor")
_COLUMNS = {"c_factor": "chl", "w_factor": "water", "element": "albino",
            "l_factor": "lignin_cellulose", "p_factor": "protein"}


def _fresnel(alpha, n0, n1):
    beta = np.arcsin(n0 / n1 * np.sin(alpha))
    plus, dif = alpha + beta, alpha - beta
    return 0.5 * (np.sin(dif) ** 2 / np.sin(plus) ** 2 + np.tan(dif) ** 2 / np.tan(plus) ** 2)


def interface_reflectances(n1, n0=1.0):
    """External (``me``) and internal (``mi``) diffuse interface reflectance
    for refractive index ``n1``, by 1-degree rectangle-rule integration, as in
    Dawson et al. (1998) / rrtm."""
    width = np.pi / 180
    alpha = np.arange(1, 91) * width
    me = 2 * np.sum(_fresnel(alpha, n0, n1) * np.sin(alpha) * np.cos(alpha) * width)
    critical = np.degrees(np.arcsin(n0 / n1))
    alpha = np.arange(1, np.floor(critical) + 1) * width  # R: seq(1, critical)
    mint = np.sum(_fresnel(alpha, n0, n1) * np.sin(alpha) * np.cos(alpha) * width)
    mi = (1 - np.sin(np.radians(critical)) ** 2) + 2 * mint
    return me, mi


def _slab_m_np(coeff):
    c = np.asarray(coeff, dtype=float)
    small = c < 1e-4
    cs = np.where(small, 1.0, c)
    full = 2 * (1 - (cs + 1) * np.exp(-cs)) / cs**2
    series = 1 - 2 * c / 3 + c**2 / 4
    return np.where(small, series, full)


def _benford_fraction(rb, tb, fraction, sqrt=np.sqrt):
    """Benford's interpolation between one and two layers (rrtm lines 117-123)."""
    top = tb ** (1 + fraction) * (((1 + tb) ** 2 - rb**2) ** (1 - fraction))
    bot1 = (1 + tb) ** (2 * (1 - fraction)) - rb**2
    bot2 = 1 + (64 / 3) * fraction * (fraction - 0.5) * (fraction - 1) * 0.001
    tif = top / (bot1 * bot2)
    rif = (1 + rb**2 - tb**2 - sqrt((1 + rb**2 - tb**2) ** 2 - 4 * rb**2 * (1 - tif**2))) / (2 * rb)
    return rif, tif


@dataclass
class Liberty(RTM):
    K: np.ndarray = field(default=None, repr=False)  # (NWL, 5), columns ordered as ABSORBERS
    me: np.ndarray = field(default=None, repr=False)
    mi: np.ndarray = field(default=None, repr=False)

    @classmethod
    def create(cls, D: float = 40.0) -> "Liberty":
        with (files("leafrtm.rtm") / "data" / "liberty.csv").open() as f:
            coef = pd.read_csv(f).iloc[:NWL]
        n1 = 1.4891 - 0.0005 * np.arange(NWL)
        me, mi = np.array([interface_reflectances(n) for n in n1]).T
        return cls(
            name="liberty",
            description="LIBERTY (Dawson et al. 1998)",
            params=("xu", "thick", "baseline", *ABSORBERS),
            wavelength=coef["wavelength"].to_numpy(dtype=float),
            fixed={"D": D},
            K=coef[[_COLUMNS[a] for a in ABSORBERS]].to_numpy(),
            me=me, mi=mi,
        )

    # -- NumPy -----------------------------------------------------------------
    def _quadratic(self, th):
        """Coefficients of a R^2 + b R + c = 0 for R_inf, plus x and Tt."""
        C = np.array([th[a] for a in ABSORBERS], dtype=float)
        coeff = th["D"] * (th["baseline"] + self.K @ C)
        me, mi, xu = self.me, self.mi, th["xu"]
        M = _slab_m_np(coeff)
        Tt = (1 - mi) * M / (1 - mi * M)
        x = xu / (1 - (1 - 2 * xu) * Tt)
        a = me * Tt + x * Tt - me - Tt - x * me * Tt
        b = 1 + x * me * Tt - 2 * x**2 * me**2 * Tt
        c = 2 * me * x**2 * Tt - x * Tt - 2 * x * me
        return a, b, c, x, Tt

    def discriminant(self, values) -> np.ndarray:
        """b^2 - 4ac per wavelength; negative where R_inf has no real root (no model)."""
        a, b, c, _, _ = self._quadratic(self.theta(values))
        return b**2 - 4 * a * c

    def validity_margin(self, values) -> float:
        """Smallest discriminant over wavelengths (weakly absorbing NIR bands go first)."""
        return float(np.min(self.discriminant(values)))

    def simulate(self, values, return_rinf=False):
        th = self.theta(values)
        me = self.me
        a, b, c, x, Tt = self._quadratic(th)
        with np.errstate(invalid="ignore"):
            R = -2 * c / (b + np.sqrt(b**2 - 4 * a * c))
            rb = 2 * x * me + x * Tt - x * Tt * 2 * x * me
            tb = np.sqrt((R - rb) * (1 - R * rb) / R)
            whole = np.floor(th["thick"])
            rif, tif = _benford_fraction(rb, tb, th["thick"] - whole)
            cur_r, cur_t = _stack_np(rb, tb, max(whole - 1, 0))
        trans = cur_t * tif / (1 - rif * cur_r)
        refl = cur_r + cur_t**2 * rif / (1 - rif * cur_r)
        return (refl, trans, R) if return_rinf else (refl, trans)

    # -- PyTensor --------------------------------------------------------------
    def simulate_pt(self, values):
        th = self.theta(values)
        C = pt.stack([pt.as_tensor_variable(th[a]) for a in ABSORBERS])
        coeff = th["D"] * (th["baseline"] + pt.dot(pt.constant(self.K), C))
        me, mi, xu = pt.constant(self.me), pt.constant(self.mi), th["xu"]
        # M(coeff) with its series near 0 (coeff >= D * baseline > 0 in practice)
        cs = pt.maximum(coeff, 1e-4)
        M = pt.switch(coeff < 1e-4, 1 - 2 * coeff / 3 + coeff**2 / 4,
                      2 * (1 - (cs + 1) * pt.exp(-cs)) / cs**2)
        Tt = (1 - mi) * M / (1 - mi * M)
        x = xu / (1 - (1 - 2 * xu) * Tt)
        a = me * Tt + x * Tt - me - Tt - x * me * Tt
        b = 1 + x * me * Tt - 2 * x**2 * me**2 * Tt
        c = 2 * me * x**2 * Tt - x * Tt - 2 * x * me
        R = -2 * c / (b + pt.sqrt(b**2 - 4 * a * c))  # NaN if no real root
        rb = 2 * x * me + x * Tt - x * Tt * 2 * x * me
        tb = pt.sqrt((R - rb) * (1 - R * rb) / R)
        whole = pt.floor(th["thick"])
        rif, tif = _benford_fraction(rb, tb, th["thick"] - whole, sqrt=pt.sqrt)
        cur_r, cur_t = _stack_pt(rb, tb, pt.maximum(whole - 1, 0))
        trans = cur_t * tif / (1 - rif * cur_r)
        refl = cur_r + cur_t**2 * rif / (1 - rif * cur_r)
        return refl, trans


def _stokes_terms(r, t, sqrt):
    D = sqrt((1 + r + t) * (1 + r - t) * (1 - r + t) * (1 - r - t))
    va = (1 + r**2 - t**2 + D) / (2 * r)
    vb = (1 - r**2 + t**2 + D) / (2 * t)
    return va, vb


def _stack_np(r, t, n):
    """(R, T) of ``n`` identical (r, t) layers, overflow-safe; n = 0 gives (0, 1)."""
    va, vb = _stokes_terms(r, t, np.sqrt)
    u = vb ** (-2 * n)
    return va * (1 - u) / (va**2 - u), vb ** (-n) * (va**2 - 1) / (va**2 - u)


def _stack_pt(r, t, n):
    va, vb = _stokes_terms(r, t, pt.sqrt)
    u = vb ** (-2 * n)
    return va * (1 - u) / (va**2 - u), vb ** (-n) * (va**2 - 1) / (va**2 - u)
