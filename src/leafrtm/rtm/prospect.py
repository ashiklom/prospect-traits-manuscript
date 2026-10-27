"""PROSPECT 4, 5, 5B and D.

A PROSPECT version is a table of specific absorption coefficients (one column
per constituent) plus a refractive index; every version shares the same
generalized plate model (``leafrtm.rtm.gpm``). The design (version = table +
constituent list, interface terms precomputed in the table) is vendored from
ashiklom/prospect-inversion-comparison @ 84504ed (``config/models.json``,
``python/prospect_np.py``). Coefficient provenance: ``data/PROVENANCE.md``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from importlib.resources import files

import numpy as np
import pandas as pd
import pytensor.tensor as pt

from leafrtm.rtm.base import RTM
from leafrtm.rtm.gpm import gpm_np, gpm_pt, interface_terms

DATA = files("leafrtm.rtm") / "data"


def _read_table(name: str) -> pd.DataFrame:
    with (DATA / f"{name}.csv").open() as f:
        return pd.read_csv(f)


# name: (table, constituents, description)
VERSIONS = {
    "prospect4": ("prospect4", ("Cab", "Cw", "Cm"), "PROSPECT-4 (Feret et al. 2008)"),
    "prospect5": ("prospect5", ("Cab", "Car", "Cw", "Cm"), "PROSPECT-5 (Feret et al. 2008)"),
    "prospect5b": ("prospect5", ("Cab", "Car", "Cbrown", "Cw", "Cm"),
                   "PROSPECT-5B (Feret et al. 2008; brown pigments, Jacquemoud et al. 2009)"),
    "prospectd": ("prospectd", ("Cab", "Car", "Canth", "Cbrown", "Cw", "Cm"),
                  "PROSPECT-D (Feret et al. 2017)"),
}


@dataclass
class Prospect(RTM):
    constituents: tuple[str, ...] = ()
    K: np.ndarray = field(default=None, repr=False)
    talf: np.ndarray = field(default=None, repr=False)
    t12: np.ndarray = field(default=None, repr=False)
    t21: np.ndarray = field(default=None, repr=False)

    @classmethod
    def from_version(cls, version: str) -> "Prospect":
        table, constituents, description = VERSIONS[version]
        coef = _read_table(table)
        talf, t12, t21 = (coef[c].to_numpy() for c in ("talf", "t12", "t21"))
        # The precomputed interface terms must agree with the refractive index.
        np.testing.assert_allclose(np.c_[talf, t12, t21],
                                   np.c_[interface_terms(coef["n"].to_numpy())], rtol=1e-8)
        return cls(name=version, description=description, params=("N", *constituents),
                   wavelength=coef["wavelength"].to_numpy(dtype=float),
                   constituents=constituents, K=coef[list(constituents)].to_numpy(),
                   talf=talf, t12=t12, t21=t21)

    def simulate(self, values, e1="fast"):
        th = self.theta(values)
        C = np.array([th[c] for c in self.constituents], dtype=float)
        k = self.K @ C / th["N"]
        return gpm_np(k, th["N"], self.talf, self.t12, self.t21, e1=e1)

    def simulate_pt(self, values, e1="fast"):
        th = self.theta(values)
        C = pt.stack([pt.as_tensor_variable(th[c]) for c in self.constituents])
        k = pt.dot(pt.constant(self.K), C) / th["N"]
        return gpm_pt(k, th["N"], pt.constant(self.talf), pt.constant(self.t12),
                      pt.constant(self.t21), e1=e1)
