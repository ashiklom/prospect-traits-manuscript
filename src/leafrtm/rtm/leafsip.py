"""LeafSIP: leaf optical properties from spectral invariants theory.

Wu, S., Zeng, Y., Hao, D., Liu, Q., Li, J., Chen, X., Asrar, G. R., Yin, G.,
Wen, J., Yang, B., Zhu, P., Chen, M. (2021). Quantifying leaf optical
properties with spectral invariants theory. Remote Sensing of Environment
253, 112131. doi:10.1016/j.rse.2020.112131

Re-implemented from the published equations and the reference MATLAB code
``SIP_Model.m`` (v1.0, 2020-11-05) at github.com/chenminbnu/SIP-RTM-Leaf
(Y. Zeng, M. Chen, S. Wu, D. Hao). As in the original, the specific absorption
coefficients are PROSPECT-D's.

    w0 = exp(-sum_i K_i C_i / (alpha * Cm))          single-scattering albedo of
                                                     a scattering element
    p  = 1 - (1 - exp(-f)) / f,   f = 2765.0 Cm      photon recollision probability
    q  = 1 - 2 exp(-g),           g = 102.8 Cm       scattering asymmetry
    w  = w0 (1 - p) / (1 - p w0)                     leaf albedo
    R  = w (1/2 + q/2 (1 - p w0) / (1 - |q| p w0))
    T  = w (1/2 - q/2 (1 - p w0) / (1 - |q| p w0))

with alpha = 600 (optimal size of the scattering element). The MATLAB code adds
``eps`` to two denominators to avoid 0/0; that is unnecessary for Cm > 0.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pytensor.tensor as pt

from leafrtm.rtm.base import RTM
from leafrtm.rtm.prospect import _read_table

ALPHA = 600.0
CONSTITUENTS = ("Cab", "Car", "Canth", "Cbrown", "Cw", "Cm")


def _sip(Kall_over_cm, Cm, exp, abs_):
    w0 = exp(-Kall_over_cm / ALPHA)
    f = 2765.0 * Cm
    g = 102.8 * Cm
    p = 1 - (1 - exp(-f)) / f
    q = 1 - 2 * exp(-g)
    w = w0 * (1 - p) / (1 - p * w0)
    ratio = (1 - p * w0) / (1 - abs_(q) * p * w0)
    return w * (0.5 + q / 2 * ratio), w * (0.5 - q / 2 * ratio)


@dataclass
class LeafSIP(RTM):
    K: np.ndarray = field(default=None, repr=False)

    @classmethod
    def create(cls) -> "LeafSIP":
        coef = _read_table("prospectd")
        return cls(name="leafsip", description="LeafSIP (Wu et al. 2021)",
                   params=CONSTITUENTS, wavelength=coef["wavelength"].to_numpy(dtype=float),
                   K=coef[list(CONSTITUENTS)].to_numpy())

    def simulate(self, values):
        th = self.theta(values)
        C = np.array([th[c] for c in CONSTITUENTS], dtype=float)
        return _sip(self.K @ C / th["Cm"], th["Cm"], np.exp, np.abs)

    def simulate_pt(self, values):
        th = self.theta(values)
        C = pt.stack([pt.as_tensor_variable(th[c]) for c in CONSTITUENTS])
        return _sip(pt.dot(pt.constant(self.K), C) / th["Cm"], th["Cm"], pt.exp, pt.abs)
