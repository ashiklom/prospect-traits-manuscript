"""Normalized schema for the harmonized leaf spectra + traits database.

One parquet file per table under ``data/processed/db/``. Hierarchy:

    dataset -> site -> plot
    entity (a physical leaf/sample; species, plot, plant)
      -> observation (an entity at one time; sample state, treatments)
           -> spectra (replicate measurements; quantity, instrument, grid)
           -> traits (leaf-level measurements)
    traits may also attach to an entity or a coarser level (``trait_level``).

Categorical metadata use the fixed enums below, so downstream analyses can
group and filter without string cleaning.
"""

from __future__ import annotations

import pandas as pd
import pandera.pandas as pa
from pandera.typing import Series

# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------
SOURCES = ["ecosis", "ecostress", "shift", "spectra_db"]
SAMPLE_STATE = ["fresh", "dried_intact", "pressed", "ground"]
TARGET_CONFIG = ["single_leaf", "leaf_stack", "needle_array", "single_needle", "unknown"]
FOREOPTIC = ["integrating_sphere", "leaf_clip", "contact_probe", "bare_fiber", "other", "unknown"]
QUANTITY = ["R", "T", "A", "PA"]
CANOPY_POSITION = ["top", "middle", "bottom"]
SUN_SHADE = ["sun", "shade"]
LEAF_AGE = ["young", "mature", "old", "senescent"]
TREATMENT_FACTOR = ["nitrogen", "water", "temperature", "co2", "light", "pathogen", "herbivory",
                    "pollutant", "genotype", "phenology_stage", "salinity", "nutrient", "other"]
GROWTH_FORM = ["tree", "shrub", "forb", "graminoid", "vine", "fern", "moss", "lichen", "aquatic",
               "succulent"]
LEAF_TYPE = ["broad", "needle", "scale", "grass"]
PHENOLOGY = ["deciduous", "evergreen"]
PS_PATHWAY = ["C3", "C4", "CAM"]
ECOSYSTEM = ["forest", "grassland", "shrubland", "savanna", "tundra", "wetland", "cropland",
             "desert", "urban", "greenhouse", "laboratory", "mixed"]
TRAIT_LEVEL = ["leaf", "plant", "plot", "species"]


def _cat(values, nullable=True):
    return pa.Field(isin=values, nullable=nullable)


class _Table(pa.DataFrameModel):
    class Config:
        coerce = True


class Datasets(_Table):
    dataset_id: Series[str] = pa.Field(unique=True)
    title: Series[str]
    source: Series[str] = _cat(SOURCES, nullable=False)
    source_id: Series[str]
    url: Series[str] = pa.Field(nullable=True)
    doi: Series[str] = pa.Field(nullable=True)
    citation: Series[str] = pa.Field(nullable=True)
    license: Series[str] = pa.Field(nullable=True)


class Sites(_Table):
    site_id: Series[str] = pa.Field(unique=True)
    dataset_id: Series[str]
    name: Series[str] = pa.Field(nullable=True)
    lat: Series[float] = pa.Field(ge=-90, le=90, nullable=True)
    lon: Series[float] = pa.Field(ge=-180, le=180, nullable=True)
    ecosystem: Series[str] = _cat(ECOSYSTEM)


class Plots(_Table):
    plot_id: Series[str] = pa.Field(unique=True)
    site_id: Series[str]
    name: Series[str] = pa.Field(nullable=True)
    lat: Series[float] = pa.Field(ge=-90, le=90, nullable=True)
    lon: Series[float] = pa.Field(ge=-180, le=180, nullable=True)


class Species(_Table):
    species_id: Series[str] = pa.Field(unique=True)
    scientific_name: Series[str]
    genus: Series[str] = pa.Field(nullable=True)
    family: Series[str] = pa.Field(nullable=True)
    growth_form: Series[str] = _cat(GROWTH_FORM)
    leaf_type: Series[str] = _cat(LEAF_TYPE)
    phenology: Series[str] = _cat(PHENOLOGY)
    ps_pathway: Series[str] = _cat(PS_PATHWAY)


class Treatments(_Table):
    treatment_id: Series[str] = pa.Field(unique=True)
    dataset_id: Series[str]
    factor: Series[str] = _cat(TREATMENT_FACTOR, nullable=False)
    level: Series[str]
    is_control: Series[pd.BooleanDtype] = pa.Field(nullable=True)


class Entities(_Table):
    entity_id: Series[str] = pa.Field(unique=True)
    dataset_id: Series[str]
    plot_id: Series[str] = pa.Field(nullable=True)
    species_id: Series[str] = pa.Field(nullable=True)
    plant_id: Series[str] = pa.Field(nullable=True)
    canopy_position: Series[str] = _cat(CANOPY_POSITION)
    sun_shade: Series[str] = _cat(SUN_SHADE)
    leaf_age: Series[str] = _cat(LEAF_AGE)


class Observations(_Table):
    observation_id: Series[str] = pa.Field(unique=True)
    entity_id: Series[str]
    dataset_id: Series[str]
    date: Series[pd.Timestamp] = pa.Field(nullable=True)
    year: Series[pd.Int64Dtype] = pa.Field(nullable=True, ge=1950, le=2030)
    sample_state: Series[str] = _cat(SAMPLE_STATE, nullable=False)
    target_config: Series[str] = _cat(TARGET_CONFIG, nullable=False)


class ObservationTreatments(_Table):
    observation_id: Series[str]
    treatment_id: Series[str]


class Instruments(_Table):
    instrument_id: Series[str] = pa.Field(unique=True)
    manufacturer: Series[str] = pa.Field(nullable=True)
    model: Series[str] = pa.Field(nullable=True)
    fwhm_nm: Series[float] = pa.Field(nullable=True, gt=0)


class Spectra(_Table):
    spectra_id: Series[str] = pa.Field(unique=True)
    observation_id: Series[str]
    quantity: Series[str] = _cat(QUANTITY, nullable=False)
    instrument_id: Series[str] = pa.Field(nullable=True)
    foreoptic: Series[str] = _cat(FOREOPTIC, nullable=False)
    grid_id: Series[str]
    qc_flags: Series[str] = pa.Field(nullable=True)
    # ``values`` (list<float32>) is validated in ``check_spectra``.


class Traits(_Table):
    observation_id: Series[str] = pa.Field(nullable=True)
    entity_id: Series[str]
    trait: Series[str]
    value: Series[float]
    unit: Series[str]
    trait_level: Series[str] = _cat(TRAIT_LEVEL, nullable=False)
    original_name: Series[str]
    original_value: Series[float] = pa.Field(nullable=True)
    original_unit: Series[str] = pa.Field(nullable=True)
    derived: Series[bool]


TABLES = {
    "datasets": Datasets, "sites": Sites, "plots": Plots, "species": Species,
    "treatments": Treatments, "entities": Entities, "observations": Observations,
    "observation_treatments": ObservationTreatments, "instruments": Instruments,
    "spectra": Spectra, "traits": Traits,
}


def check_spectra(spectra: pd.DataFrame, grids: pd.DataFrame) -> None:
    """Every spectrum's values match its grid length; grids strictly increase."""
    lengths = grids.set_index("grid_id")["wavelength"].map(len)
    bad = spectra["values"].map(len).to_numpy() != lengths.reindex(spectra["grid_id"]).to_numpy()
    if bad.any():
        raise ValueError(f"{bad.sum()} spectra do not match their grid length")
    for gid, wl in zip(grids["grid_id"], grids["wavelength"]):
        if len(wl) > 1 and not (pd.Series(wl).diff().iloc[1:] > 0).all():
            raise ValueError(f"Grid {gid} is not strictly increasing")


def check_links(db: dict[str, pd.DataFrame]) -> dict[str, int]:
    """Count orphan foreign keys (should all be zero)."""
    def orphans(child, key, parent, pkey=None):
        c = db[child][key].dropna()
        return int((~c.isin(db[parent][pkey or key])).sum())

    return {
        "sites.dataset_id": orphans("sites", "dataset_id", "datasets"),
        "plots.site_id": orphans("plots", "site_id", "sites"),
        "entities.plot_id": orphans("entities", "plot_id", "plots"),
        "entities.species_id": orphans("entities", "species_id", "species"),
        "observations.entity_id": orphans("observations", "entity_id", "entities"),
        "spectra.observation_id": orphans("spectra", "observation_id", "observations"),
        "spectra.grid_id": orphans("spectra", "grid_id", "wavelength_grids"),
        "traits.entity_id": orphans("traits", "entity_id", "entities"),
        "traits.observation_id": orphans("traits", "observation_id", "observations"),
        "observation_treatments.treatment_id": orphans("observation_treatments", "treatment_id",
                                                       "treatments"),
    }


def validate(db: dict[str, pd.DataFrame]) -> None:
    for name, model in TABLES.items():
        model.validate(db[name], lazy=True)
    check_spectra(db["spectra"], db["wavelength_grids"])
    bad = {k: v for k, v in check_links(db).items() if v}
    if bad:
        raise ValueError(f"Orphan keys: {bad}")
