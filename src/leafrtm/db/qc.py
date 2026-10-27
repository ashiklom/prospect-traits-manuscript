"""Quality-control flags. Flags are recorded; values are never modified.

Spectral flags (``spectra.qc_flags``, ';'-separated):

* ``out_of_range``: more than 5 % of 400-2500 nm bands outside [-0.02, 1.02]
  (R, T, A)
* ``flat``: standard deviation across 400-2500 nm below 0.002
* ``sparse``: fewer than 50 bands in 400-2500 nm
* ``jump``: a step at a detector join (~1000 or ~1800 nm) larger than 3x the
  local band-to-band variability and > 0.01
* ``replicate_outlier``: RMS difference from the observation's replicate
  median above 4 x the dataset's median RMS (and > 0.01)

Trait screening: values outside physically plausible bounds are dropped
(listed in the harmonization report).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

PLAUSIBLE = {  # canonical units (vocab/traits.yaml)
    "chl_area": (0, 150), "chla_area": (0, 120), "chlb_area": (0, 60), "car_area": (0, 50),
    "anth_area": (0, 100), "lma": (2, 1500), "ewt": (1, 2000), "lwc": (1, 99), "ldmc": (20, 900),
    "n_mass": (0.1, 10), "n_area": (0.01, 20), "c_mass": (20, 65), "cn_ratio": (3, 300),
    "p_mass": (0.005, 2), "lignin_mass": (0, 60), "cellulose_mass": (0, 80),
    "vcmax_area": (1, 500), "vcmax25_area": (1, 500), "jmax_area": (1, 800),
    "jmax25_area": (1, 800), "leaf_thickness": (20, 5000), "leaf_area": (0.001, 10000),
}


def screen_traits(traits: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Drop implausible values; returns (kept, dropped)."""
    lo = traits["trait"].map(lambda t: PLAUSIBLE.get(t, (-np.inf, np.inf))[0])
    hi = traits["trait"].map(lambda t: PLAUSIBLE.get(t, (-np.inf, np.inf))[1])
    ok = (traits["value"] >= lo) & (traits["value"] <= hi) & np.isfinite(traits["value"])
    return traits[ok].copy(), traits[~ok].copy()


def spectral_flags(values: np.ndarray, wavelength: np.ndarray, quantity: str) -> list[str]:
    flags = []
    m = (wavelength >= 400) & (wavelength <= 2500)
    v = values[m]
    w = wavelength[m]
    if len(v) < 50:
        return ["sparse"]
    if quantity in ("R", "T", "A"):
        if np.mean((v < -0.02) | (v > 1.02)) > 0.05:
            flags.append("out_of_range")
    if np.std(v) < 0.002:
        flags.append("flat")
    for join in (1000.0, 1800.0):
        i = np.searchsorted(w, join)
        if 5 <= i < len(w) - 5 and w[i] - w[i - 1] <= 3:
            local = np.median(np.abs(np.diff(v[i - 5:i + 5])))
            step = abs(v[i] - v[i - 1])
            if step > 0.01 and step > 3 * max(local, 1e-4):
                flags.append("jump")
                break
    return flags


def replicate_outliers(spectra: pd.DataFrame, grids: dict) -> pd.Series:
    """Boolean Series (index of ``spectra``): replicate far from its group median."""
    rms = pd.Series(np.nan, index=spectra.index)
    for _, g in spectra.groupby(["observation_id", "quantity", "grid_id"]):
        if len(g) < 3:
            continue
        X = np.stack(g["values"].to_numpy()).astype(float)
        med = np.median(X, axis=0)
        rms.loc[g.index] = np.sqrt(np.mean((X - med) ** 2, axis=1))
    dataset = spectra["observation_id"].str.split(":").str[0]
    typical = rms.groupby(dataset).transform("median")
    return ((rms > 4 * typical) & (rms > 0.01)).fillna(False).astype(bool)
