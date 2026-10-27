"""Select observations to fit and build their inversion inputs.

Each observation becomes a list of ``obsop.Block`` objects (replicate means
of one quantity on one grid). Replicates flagged as QC outliers are left out
of the mean; spectra flagged out_of_range/flat/sparse are not used.
Observed bands are never resampled.
"""

from __future__ import annotations

import pickle
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

from leafrtm.db.harmonize import load
from leafrtm.obsop import Block

EXCLUDE_FLAGS = ("out_of_range", "flat", "sparse", "replicate_outlier")
SRF_STEP = 8.0  # nm: band-averaged data coarser than this is matched with a Gaussian SRF


def usable_spectra(spectra: pd.DataFrame) -> pd.DataFrame:
    f = spectra["qc_flags"].fillna("")
    bad = np.zeros(len(spectra), bool)
    for flag in EXCLUDE_FLAGS:
        bad |= f.str.contains(flag).to_numpy()
    return spectra[~bad]


def observation_table(db=None, quantities=("R", "T", "A", "PA")) -> pd.DataFrame:
    """One row per observation with fit-relevant attributes."""
    db = db or load()
    s = usable_spectra(db["spectra"])
    s = s[s["quantity"].isin(quantities)]
    q = s.groupby("observation_id")["quantity"].agg(lambda x: "".join(sorted(set(x))))
    obs = db["observations"].set_index("observation_id").join(q.rename("quantities"), how="inner")
    obs = obs.join(db["entities"].set_index("entity_id")[["species_id", "plant_id"]], on="entity_id")
    tr = db["traits"]
    n_traits = tr.dropna(subset=["observation_id"]).groupby("observation_id")["trait"].nunique()
    obs["n_traits"] = n_traits.reindex(obs.index).fillna(0).astype(int)
    obs["has_RT"] = obs["quantities"].str.contains("R") & obs["quantities"].str.contains("T")
    return obs.reset_index()


def batch_of(observation_id: str, n_batches: int) -> int:
    return zlib.crc32(observation_id.encode()) % n_batches


def select(obs: pd.DataFrame, fitset: dict, seed=20261027) -> pd.DataFrame:
    """Fit set: all observations (full) or a stratified pilot subset."""
    obs = obs[obs["quantities"].str.contains("|".join(fitset.get("quantities", ["R", "T"])))]
    if fitset.get("datasets"):
        obs = obs[obs["dataset_id"].isin(fitset["datasets"])]
    n = fitset.get("n_observations")
    if n:
        rng = np.random.default_rng(seed)
        strata = fitset.get("strata", ["dataset_id", "sample_state", "has_RT"])

        def draw(pool, k):
            if k <= 0 or pool.empty:
                return pool.iloc[:0]
            w = pool.groupby(strata)["observation_id"].transform("size")
            # equal weight per stratum, favouring observations with measured traits
            p = (1.0 / w) * (1 + np.minimum(pool["n_traits"], 5))
            idx = rng.choice(len(pool), size=min(k, len(pool)), replace=False, p=p / p.sum())
            return pool.iloc[idx]

        # R+T leaves carry the R-only vs R+T and held-out-T comparisons: guarantee a quota
        rt = draw(obs[obs["has_RT"]], fitset.get("min_rt", 0))
        rest = draw(obs[~obs["observation_id"].isin(rt["observation_id"])], n - len(rt))
        obs = obs[obs["observation_id"].isin(pd.concat([rt, rest])["observation_id"])]
    obs = obs.copy()
    obs["batch"] = [batch_of(o, fitset["n_batches"]) for o in obs["observation_id"]]
    return obs


def blocks_for(obs_ids, db=None, mask_ranges=(), quantities=("R", "T", "A", "PA")) -> dict:
    """{observation_id: [Block, ...]} with replicate means per quantity and grid."""
    db = db or load()
    grids = {g: (np.asarray(w, float), np.asarray(s, bool)) for g, w, s in zip(
        db["wavelength_grids"]["grid_id"], db["wavelength_grids"]["wavelength"],
        db["wavelength_grids"]["segment_start"])}
    spec = usable_spectra(db["spectra"])
    spec = spec[spec["observation_id"].isin(set(obs_ids)) & spec["quantity"].isin(quantities)]
    spec = spec.merge(db["instruments"][["instrument_id"]], on="instrument_id", how="left")
    out = {}
    # replicates are averaged only within one instrument and grid (some
    # observations combine spectrometers with different radiometric offsets)
    spec = spec.assign(instrument_id=spec["instrument_id"].fillna("unknown"))
    for (oid, q, gid, _), g in spec.groupby(["observation_id", "quantity", "grid_id",
                                              "instrument_id"]):
        wl, seg = grids[gid]
        X = np.stack(g["values"].to_numpy()).astype(float)
        keep = np.ones(len(wl), bool)
        for lo, hi in mask_ranges:
            keep &= ~((wl >= lo) & (wl <= hi))
        seg = seg.copy()
        # a band following a masked gap starts a new segment only at detector joins;
        # gaps decorrelate through distance
        step = np.median(np.diff(wl)) if len(wl) > 1 else 1.0
        y = X.mean(axis=0)
        out.setdefault(oid, []).append(Block(
            quantity=q, wavelength=wl[keep], y=y[keep], n_rep=len(X),
            rep_sd=X.std(axis=0, ddof=1)[keep] if len(X) > 1 else None,
            segment_start=seg[keep], fwhm=step if step > SRF_STEP else None,
            instrument_class=str(g["foreoptic"].iloc[0])))
    return out


def data_use_variants(blocks: list[Block], variants=("all",)) -> dict[str, list[Block]]:
    """Subsets of an observation's blocks by quantity.

    ``all`` uses every block; otherwise a variant names quantities joined by
    '+' (e.g. ``R``, ``R+T``, ``A``, ``PA``) and exists only if the
    observation has all of them. ``R`` additionally requires another
    quantity to exist (so it can be predicted out of sample).
    """
    out = {}
    qs = {b.quantity for b in blocks}
    for v in variants:
        if v == "all":
            out[v] = blocks
            continue
        want = set(v.split("+"))
        if not want <= qs or (v == "R" and len(qs) == 1):
            continue
        out[v] = [b for b in blocks if b.quantity in want]
    return out


def write_batch(path: Path, blocks: dict, meta: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump({"blocks": blocks, "meta": meta}, f, protocol=pickle.HIGHEST_PROTOCOL)


def read_batch(path: Path) -> dict:
    with open(path, "rb") as f:
        return pickle.load(f)
