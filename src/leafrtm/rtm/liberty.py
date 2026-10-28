"""LIBERTY (Leaf Incorporating Biochemistry Exhibiting Reflectance and
Transmittance Yields; Dawson, Curran & Plummer 1998, RSE 65:50-60,
doi:10.1016/S0034-4257(98)00007-8), in two variants:

* ``liberty`` (class ``Liberty``): fixes a bug in the original paper. Its
  infinite-reflectance equation (Eq. 16) is not the sum of Eq. 15 it is derived
  from; ``liberty`` uses the correct sum. Use this one.
* ``liberty_dawson`` (class ``LibertyDawson``): the model as published by Dawson
  et al. (1998) and implemented in rrtm, kept for provenance and comparison. It
  has no real solution for weakly absorbing, strongly scattering leaves.

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

One correction relative to rrtm, in both variants: the internal reflectance
``mi`` (Dawson et al. Eq. 2) integrates the Fresnel reflectance of light inside
the cell at internal angles up to the critical angle. rrtm evaluates the
reflectance of light arriving from *outside* at those angles (it pairs each
angle with its refraction into the cell, not out of it), and sums whole degrees
only, missing the rise to total internal reflection just below the critical
angle; at 800 nm it gives mi = 0.541 instead of 0.564. Eqs. 1 and 2 have closed
forms in Stern's (1964) diffuse transmissivity ``tav(90, n)``: ``me = 1 - tav``
and ``mi = 1 - tav / n^2``, which we use (rrtm's 1-degree sum for ``me`` was
within 6e-5 of it). Eq. 1 as printed lacks the factor 2 that normalizes the
cosine-weighted average; Eq. 2 and rrtm both have it.

The bug fixed in ``liberty``. Eq. 15 sums the radiation leaving a surface layer of cells
lying on leaf material of reflectance R (Fig. 1). Its sum is

    R = 2 x me + x (1 - 2 x me) tau (1 - me R) / [(1 - me R) - (1 - x)(1 - me) tau R],

but Eq. 16 places ``2 x me`` over the common denominator without multiplying it
by that denominator. rrtm (and ``liberty_dawson``) implement Eq. 16. Eq. 15 is the
adding equation ``R = r_d + t_d t_u R / (1 - r_u R)`` for the surface layer of
Fig. 1, with

    r_d = 2 x me + x (1 - 2 x me) tau     (Eq. 17)    t_d = (1 - x)(1 - 2 x me) tau
    r_u = me + (1 - x)(1 - me) tau                    t_u = x (1 - me) tau

(d: light arriving from above; u: from below). The layer conserves energy from
either side (it absorbs (1 - 2 x me)(1 - tau) and (1 - me)(1 - tau)), so
sqrt(r_d r_u) + sqrt(t_d t_u) <= 1 and R_inf is real for any absorption, tending
to 1 as absorption vanishes. Eq. 16 instead corresponds to a layer that returns
more light than it receives when absorption is weak (me (1 + me) / 2 more at
zero absorption), which is why it has no real root there. ``liberty``
replaces Eq. 16 with the sum of Eq. 15 and changes nothing else (Eqs. 1-9, 17,
19 and Benford's thickness interpolation). It evaluates Eq. 19 as
``T^2 = t_d t_u (1 - R r_d) / (1 - R r_u)``, equal to ``(R - r)(1 - R r) / R`` by
the adding equation but free of the cancellation in ``R - r`` under strong
absorption.

Notes for inference: ``D`` (cell diameter) multiplies every absorption term,
so it is exactly non-identifiable from the concentrations; it is held fixed
(``fixed``) at Dawson's typical value. ``thick`` must be >= 1 (the model jumps
from 2 to 1 layers at thick = 1). In ``liberty_dawson``, large ``xu`` with
weak absorption leaves the quadratic without a real root (rrtm diverges); the
PyTensor version returns NaN there, which the sampler rejects. ``liberty`` is
defined for all parameter values.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from importlib.resources import files

import numpy as np
import pandas as pd
import pytensor.tensor as pt

from leafrtm.rtm.base import RTM
from leafrtm.rtm.gpm import tav

NWL = 420
ABSORBERS = ("c_factor", "w_factor", "element", "l_factor", "p_factor")
_COLUMNS = {"c_factor": "chl", "w_factor": "water", "element": "albino",
            "l_factor": "lignin_cellulose", "p_factor": "protein"}


def interface_reflectances(n):
    """Diffuse reflectance of the cell surface for refractive index ``n``:
    ``me`` for light arriving from the air (Eq. 1) and ``mi`` for light inside
    the cell (Eq. 2)."""
    t12 = tav(90, n)
    return 1 - t12, 1 - t12 / np.asarray(n, dtype=float) ** 2


def _slab_m_np(coeff):
    c = np.asarray(coeff, dtype=float)
    small = c < 1e-4
    cs = np.where(small, 1.0, c)
    full = 2 * (1 - (cs + 1) * np.exp(-cs)) / cs**2
    series = 1 - 2 * c / 3 + c**2 / 4
    return np.where(small, series, full)


def _slab_m_pt(coeff):
    # M(coeff) with its series near 0 (coeff >= D * baseline > 0 in practice)
    cs = pt.maximum(coeff, 1e-4)
    return pt.switch(coeff < 1e-4, 1 - 2 * coeff / 3 + coeff**2 / 4,
                     2 * (1 - (cs + 1) * pt.exp(-cs)) / cs**2)


def _cell(mi, xu, M):
    """Cell transmission tau (Eq. 6) and upward-scattering coefficient x (Eq. 9)."""
    Tt = (1 - mi) * M / (1 - mi * M)
    return xu / (1 - (1 - 2 * xu) * Tt), Tt


def _surface_layer(me, x, Tt):
    """Reflectance and transmittance of the surface cell layer of Fig. 1 for
    light from above (r_d, t_d) and from below (r_u, t_u)."""
    rd = 2 * x * me + x * Tt - x * Tt * 2 * x * me
    td = (1 - x) * (1 - 2 * x * me) * Tt
    ru = me + (1 - x) * (1 - me) * Tt
    tu = x * (1 - me) * Tt
    return rd, td, ru, tu


def _benford_fraction(rb, tb, fraction, sqrt=np.sqrt):
    """Benford's interpolation between one and two layers (rrtm lines 117-123)."""
    top = tb ** (1 + fraction) * (((1 + tb) ** 2 - rb**2) ** (1 - fraction))
    bot1 = (1 + tb) ** (2 * (1 - fraction)) - rb**2
    bot2 = 1 + (64 / 3) * fraction * (fraction - 0.5) * (fraction - 1) * 0.001
    tif = top / (bot1 * bot2)
    rif = (1 + rb**2 - tb**2 - sqrt((1 + rb**2 - tb**2) ** 2 - 4 * rb**2 * (1 - tif**2))) / (2 * rb)
    return rif, tif


def _leaf(rb, tb, thick, xp):
    """(R, T) of a leaf ``thick`` layers deep; ``xp`` is ``np`` or ``pt``."""
    whole = xp.floor(thick)
    rif, tif = _benford_fraction(rb, tb, thick - whole, sqrt=xp.sqrt)
    cur_r, cur_t = _stack(rb, tb, xp.maximum(whole - 1, 0), xp.sqrt)
    trans = cur_t * tif / (1 - rif * cur_r)
    refl = cur_r + cur_t**2 * rif / (1 - rif * cur_r)
    return refl, trans


@dataclass
class LibertyDawson(RTM):
    """LIBERTY as published by Dawson et al. (1998), Eq. 16 included (as in rrtm)."""

    K: np.ndarray = field(default=None, repr=False)  # (NWL, 5), columns ordered as ABSORBERS
    me: np.ndarray = field(default=None, repr=False)
    mi: np.ndarray = field(default=None, repr=False)

    _label = ("liberty_dawson", "LIBERTY as published (Dawson et al. 1998)")

    @classmethod
    def create(cls, D: float = 40.0) -> "LibertyDawson":
        with (files("leafrtm.rtm") / "data" / "liberty.csv").open() as f:
            coef = pd.read_csv(f).iloc[:NWL]
        me, mi = interface_reflectances(1.4891 - 0.0005 * np.arange(NWL))
        name, description = cls._label
        return cls(
            name=name,
            description=description,
            params=("xu", "thick", "baseline", *ABSORBERS),
            wavelength=coef["wavelength"].to_numpy(dtype=float),
            fixed={"D": D},
            K=coef[[_COLUMNS[a] for a in ABSORBERS]].to_numpy(),
            me=me, mi=mi,
        )

    def _coefficients(self, me, x, Tt):
        """a, b, c of a R^2 + b R + c = 0 for R_inf (Eq. 16)."""
        a = me * Tt + x * Tt - me - Tt - x * me * Tt
        b = 1 + x * me * Tt - 2 * x**2 * me**2 * Tt
        c = 2 * me * x**2 * Tt - x * Tt - 2 * x * me
        return a, b, c

    def _layer(self, me, x, Tt, xp):
        """Single-layer reflectance (Eq. 17), transmittance (Eq. 19) and R_inf."""
        a, b, c = self._coefficients(me, x, Tt)
        R = -2 * c / (b + xp.sqrt(b**2 - 4 * a * c))  # NaN if no real root
        rb = _surface_layer(me, x, Tt)[0]
        tb = xp.sqrt((R - rb) * (1 - R * rb) / R)
        return rb, tb, R

    # -- NumPy -----------------------------------------------------------------
    def _cell_np(self, th):
        C = np.array([th[a] for a in ABSORBERS], dtype=float)
        coeff = th["D"] * (th["baseline"] + self.K @ C)
        return _cell(self.mi, th["xu"], _slab_m_np(coeff))

    def discriminant(self, values) -> np.ndarray:
        """b^2 - 4ac per wavelength; negative where R_inf has no real root (no model)."""
        a, b, c = self._coefficients(self.me, *self._cell_np(self.theta(values)))
        return b**2 - 4 * a * c

    def validity_margin(self, values) -> float:
        """Smallest discriminant over wavelengths (weakly absorbing NIR bands go first)."""
        return float(np.min(self.discriminant(values)))

    def simulate(self, values, return_rinf=False):
        th = self.theta(values)
        x, Tt = self._cell_np(th)
        with np.errstate(invalid="ignore"):
            rb, tb, R = self._layer(self.me, x, Tt, np)
            refl, trans = _leaf(rb, tb, th["thick"], np)
        return (refl, trans, R) if return_rinf else (refl, trans)

    # -- PyTensor --------------------------------------------------------------
    def simulate_pt(self, values):
        th = self.theta(values)
        C = pt.stack([pt.as_tensor_variable(th[a]) for a in ABSORBERS])
        coeff = th["D"] * (th["baseline"] + pt.dot(pt.constant(self.K), C))
        x, Tt = _cell(pt.constant(self.mi), th["xu"], _slab_m_pt(coeff))
        rb, tb, _ = self._layer(pt.constant(self.me), x, Tt, pt)
        return _leaf(rb, tb, th["thick"], pt)


@dataclass
class Liberty(LibertyDawson):
    """LIBERTY with a bug in the original paper fixed.

    Dawson et al. (1998) Eq. 16, the infinite-reflectance equation, is not the
    sum of Eq. 15 it is derived from: ``2 x me`` is put over the common
    denominator without being multiplied by it. This class uses the correct sum
    (the adding equation of an energy-conserving surface layer), so R_inf is
    real for all parameter values; everything else is as in ``LibertyDawson``.
    See the module docstring for the derivation."""

    _label = ("liberty", "LIBERTY with Eq. 16 corrected (Dawson et al. 1998 bug fixed)")

    def _coefficients(self, me, x, Tt):
        """a, b, c of a R^2 + b R + c = 0 for R_inf (Eq. 15 summed)."""
        rd, td, ru, tu = _surface_layer(me, x, Tt)
        return -ru, 1 + rd * ru - td * tu, -rd

    def _layer(self, me, x, Tt, xp):
        rd, td, ru, tu = _surface_layer(me, x, Tt)
        b = 1 + rd * ru - td * tu
        # b^2 - 4 ru rd >= 0 analytically, = 0 only without absorption; clip rounding
        R = 2 * rd / (b + xp.sqrt(xp.maximum(b**2 - 4 * ru * rd, 0)))
        tb = xp.sqrt(td * tu * (1 - R * rd) / (1 - R * ru))
        return rd, tb, R

    validity_margin = RTM.validity_margin  # defined for all parameter values


def _stokes_terms(r, t, sqrt):
    D = sqrt((1 + r + t) * (1 + r - t) * (1 - r + t) * (1 - r - t))
    va = (1 + r**2 - t**2 + D) / (2 * r)
    vb = (1 - r**2 + t**2 + D) / (2 * t)
    return va, vb


def _stack(r, t, n, sqrt):
    """(R, T) of ``n`` identical (r, t) layers, overflow-safe; n = 0 gives (0, 1)."""
    va, vb = _stokes_terms(r, t, sqrt)
    u = vb ** (-2 * n)
    return va * (1 - u) / (va**2 - u), vb ** (-n) * (va**2 - 1) / (va**2 - u)
