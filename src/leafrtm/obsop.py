"""Observation operator: maps RTM output on its native grid onto observed
spectra, without resampling the observations.

An observation (one leaf at one time) has one or more *blocks*. A block is
the replicate mean of all spectra of one quantity measured on one wavelength
grid. Quantities:

* ``R``  reflectance, ``T`` transmittance (directly modeled)
* ``A``  absorbance ``1 - R - T``
* ``PA`` pseudo-absorbance ``log10(1/R)``

For each observed band the operator stores a few (index, weight) pairs into
the model's native grid. Depending on the instrument this is index selection
(observed band on the model grid), linear interpolation, or a Gaussian
spectral response function when the instrument's FWHM is clearly coarser than
the RTM's calibration data.

All blocks of an observation are concatenated into flat arrays (see
``ObsData``) so the likelihood is a handful of vectorized gathers.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

QUANTITIES = ("R", "T", "A", "PA")
QCODE = {q: i for i, q in enumerate(QUANTITIES)}
# Instrument FWHM (nm) above which the model is convolved with a Gaussian SRF
# instead of interpolated. RTM coefficients were calibrated on spectrometer
# data with a few nm of smoothing, so convolving finer instruments would
# double-smooth.
SRF_FWHM_THRESHOLD = 8.0


@dataclass
class Block:
    """Replicate-mean spectrum of one quantity on one wavelength grid."""

    quantity: str
    wavelength: np.ndarray
    y: np.ndarray
    n_rep: int = 1
    rep_sd: np.ndarray | None = None
    segment_start: np.ndarray | None = None  # True at detector/segment boundaries
    fwhm: float | None = None
    instrument_class: str = "unknown"

    def __post_init__(self):
        self.wavelength = np.asarray(self.wavelength, dtype=float)
        self.y = np.asarray(self.y, dtype=float)
        if self.quantity not in QCODE:
            raise ValueError(f"Unknown quantity {self.quantity!r}")
        if self.segment_start is None:
            self.segment_start = np.zeros(len(self.wavelength), dtype=bool)
        self.segment_start = np.asarray(self.segment_start, dtype=bool).copy()
        if len(self.wavelength):
            self.segment_start[0] = True
        if np.any(np.diff(self.wavelength) <= 0):
            raise ValueError("Block wavelengths must be strictly increasing")


def clean_grid(wavelength, values, segment_start=None):
    """Sort bands, average exact duplicates, drop non-finite values.

    ``values`` may be 1-D or (n_bands, n_replicates). Returns (wavelength,
    values, segment_start) with strictly increasing wavelengths.
    """
    wl = np.asarray(wavelength, dtype=float)
    v = np.asarray(values, dtype=float)
    v2 = v.reshape(len(wl), -1)
    seg = np.zeros(len(wl), bool) if segment_start is None else np.asarray(segment_start, bool)
    keep = np.isfinite(wl) & np.all(np.isfinite(v2), axis=1)
    wl, v2, seg = wl[keep], v2[keep], seg[keep]
    order = np.argsort(wl, kind="stable")
    wl, v2, seg = wl[order], v2[order], seg[order]
    uniq, inv = np.unique(wl, return_inverse=True)
    if len(uniq) < len(wl):
        sums = np.zeros((len(uniq), v2.shape[1]))
        np.add.at(sums, inv, v2)
        v2 = sums / np.bincount(inv)[:, None]
        seg = np.bincount(inv, weights=seg)[...] > 0
        wl = uniq
    return wl, (v2[:, 0] if v.ndim == 1 else v2), seg


def operator(model_wl, obs_wl, fwhm=None, threshold=SRF_FWHM_THRESHOLD):
    """(idx, w) arrays of shape (n_obs, k) such that ``pred = sum(w * m[idx], 1)``.

    Observed bands outside the model grid get weight 0 (callers should mask
    them; ``make_obsdata`` drops them).
    """
    model_wl = np.asarray(model_wl, dtype=float)
    obs_wl = np.asarray(obs_wl, dtype=float)
    if fwhm is not None and fwhm > threshold:
        sd = fwhm / (2 * np.sqrt(2 * np.log(2)))
        half = int(np.ceil(3 * sd / np.min(np.diff(model_wl))))
        centre = np.searchsorted(model_wl, obs_wl)
        offs = np.arange(-half, half + 1)
        idx = np.clip(centre[:, None] + offs[None, :], 0, len(model_wl) - 1)
        w = np.exp(-0.5 * ((model_wl[idx] - obs_wl[:, None]) / sd) ** 2)
        # duplicate clipped indices must not be double counted
        dup = np.zeros_like(w, dtype=bool)
        dup[:, 1:] = idx[:, 1:] == idx[:, :-1]
        w[dup] = 0
        w /= w.sum(axis=1, keepdims=True)
        return idx, w
    j = np.clip(np.searchsorted(model_wl, obs_wl, side="right") - 1, 0, len(model_wl) - 2)
    x0, x1 = model_wl[j], model_wl[j + 1]
    frac = np.clip((obs_wl - x0) / (x1 - x0), 0.0, 1.0)
    idx = np.c_[j, j + 1]
    w = np.c_[1 - frac, frac]
    inside = (obs_wl >= model_wl[0]) & (obs_wl <= model_wl[-1])
    w[~inside] = 0
    return idx, w


def model_quantities_np(R, T):
    """Stack modeled R, T, A, PA (rows ordered as ``QUANTITIES``)."""
    with np.errstate(divide="ignore"):
        return np.stack([R, T, 1 - R - T, -np.log10(R)])


@dataclass
class ObsData:
    """Flat, concatenated representation of all blocks of one observation.

    Arrays have one entry per observed band (``n``), except ``idx``/``w``
    (n, k). ``flat_idx`` indexes the flattened (4, n_native) quantity stack.
    """

    y: np.ndarray
    qcode: np.ndarray
    block: np.ndarray
    flat_idx: np.ndarray
    w: np.ndarray
    dlam: np.ndarray  # distance to previous band in the same block (nm); 0 at starts
    start: np.ndarray  # True at block/segment starts
    wavelength: np.ndarray
    pair: np.ndarray  # for T bands: index of the R band at the same wavelength, else -1
    n_native: int
    blocks: list = field(default_factory=list, repr=False)

    @property
    def quantities(self) -> tuple[str, ...]:
        return tuple(q for q in QUANTITIES if np.any(self.qcode == QCODE[q]))

    def __len__(self):
        return len(self.y)

    def predict_np(self, R, T):
        Q = model_quantities_np(R, T).ravel()
        return np.sum(self.w * Q[self.flat_idx], axis=1)


def make_obsdata(blocks, model_wl) -> ObsData:
    """Concatenate blocks, dropping bands outside the model's wavelength range.

    Within each block, observed bands are kept as measured (no resampling).
    """
    model_wl = np.asarray(model_wl, dtype=float)
    lo, hi = model_wl[0], model_wl[-1]
    parts = []
    kmax = 0
    kept = []
    for b in blocks:
        inside = (b.wavelength >= lo) & (b.wavelength <= hi)
        if not inside.any():
            continue
        wl = b.wavelength[inside]
        start = b.segment_start[inside].copy()
        start[0] = True
        # a band whose predecessor was dropped starts a new segment
        dropped_before = np.r_[True, ~inside[:-1]][inside]
        start |= dropped_before
        idx, w = operator(model_wl, wl, b.fwhm)
        kmax = max(kmax, idx.shape[1])
        parts.append((b, wl, b.y[inside], start, idx, w))
        kept.append(b)
    if not parts:
        raise ValueError("No observed bands inside the model's wavelength range")

    def pad(a, k, fill):
        return np.pad(a, ((0, 0), (0, k - a.shape[1])), constant_values=fill)

    y, qc, blk, fidx, ww, dl, st, wls = [], [], [], [], [], [], [], []
    n_native = len(model_wl)
    for bi, (b, wl, yb, start, idx, w) in enumerate(parts):
        q = QCODE[b.quantity]
        y.append(yb)
        qc.append(np.full(len(wl), q))
        blk.append(np.full(len(wl), bi))
        fidx.append(pad(idx, kmax, 0) + q * n_native)
        ww.append(pad(w, kmax, 0.0))
        d = np.r_[0.0, np.diff(wl)]
        d[start] = 0.0
        dl.append(d)
        st.append(start)
        wls.append(wl)
    out = ObsData(
        y=np.concatenate(y), qcode=np.concatenate(qc), block=np.concatenate(blk),
        flat_idx=np.concatenate(fidx).astype(np.int64), w=np.concatenate(ww),
        dlam=np.concatenate(dl), start=np.concatenate(st), wavelength=np.concatenate(wls),
        pair=np.full(sum(len(p[1]) for p in parts), -1, dtype=np.int64),
        n_native=n_native, blocks=kept,
    )
    # pair T bands with R bands on identical wavelengths (for correlated R/T errors)
    r_pos = {(round(wl, 6)): i for i, (wl, q) in enumerate(zip(out.wavelength, out.qcode))
             if q == QCODE["R"]}
    for i, (wl, q) in enumerate(zip(out.wavelength, out.qcode)):
        if q == QCODE["T"]:
            out.pair[i] = r_pos.get(round(wl, 6), -1)
    # the separable R/T model needs identical grids: unpair partially paired blocks
    for bi in np.unique(out.block[out.qcode == QCODE["T"]]):
        sel = out.block == bi
        if np.any(out.pair[sel] < 0):
            out.pair[sel] = -1
    return out
