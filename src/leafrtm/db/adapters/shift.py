"""Load NASA SHIFT dried & ground leaf spectra with their sample traits.

Spectra (6 scans per sample and preparation) come from ORNL DAAC 2244;
fresh-leaf LMA, LWC and chlorophyll (CCM-300) from 2233; elemental and
fiber chemistry of the oven-dried material from 2337; plot and species
metadata from 2295. Everything is joined on ``sample_number``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from leafrtm.db.adapters.core import DatasetTables, to_tables
from leafrtm.db.sources import shift
from leafrtm.db.sources.common import raw_dir


def _read(ds, name):
    p = raw_dir("shift") / ds / name
    if not p.exists():
        shift.fetch_all()
    return pd.read_csv(p, na_values=["NA", -9999, "-9999"], low_memory=False)


def load(cfg: dict) -> DatasetTables:
    spectra = []
    for prep, f in (("oven_dried", "Ovendried_Spectra.csv"), ("freeze_dried", "Flashfrozen_Spectra.csv")):
        df = _read("SHIFT_DriedGround_Leaf_Reflec", f)
        spectra.append(df.assign(prep=prep))
    spec = pd.concat(spectra, ignore_index=True)
    spec["sample_number"] = pd.to_numeric(spec["sample_number"], errors="coerce")
    wl_cols = [c for c in spec.columns if c.isdigit()]
    meta = spec.drop(columns=wl_cols)
    vals = spec[wl_cols].to_numpy(float)
    wl = np.array([float(c) for c in wl_cols])

    traits = _read("SHIFT_Leaf_Traits_Chl_SB_CA", "SHIFT_Leaf_Traits_LMA_LWC_Chl.csv")
    elem = _read("SHIFT_Foliar_Chemical_Analysis", "OD_elemental_daac.csv")
    fiber = _read("SHIFT_Foliar_Chemical_Analysis", "OD_AnkomFiber_daac.csv")
    samples = _read("SHIFT_Field_Survey_Metadata", "SHIFT_vegetation_sample_list.csv")
    samples.columns = [c.lower() for c in samples.columns]
    samples["sample_number"] = pd.to_numeric(samples["sample_number"], errors="coerce")

    info = samples[["sample_number", "plot_event_id", "species_or_type", "sample_date",
                    "plot_name", "phenophase"]].drop_duplicates("sample_number")
    tcols = ["sample_number", "latitude", "longitude", "LWC", "CFR", "CHL", "LMA"]
    ecols = ["sample_number", "Nitrogen_perc", "Phosphorus_perc", "Magnesium_perc",
             "Potassium_perc", "Calcium_perc", "Sulfur_perc"]
    fcols = ["sample_number", "cellulose_perc", "lignin_perc"]
    for d in (traits, elem, fiber):
        d["sample_number"] = pd.to_numeric(d["sample_number"], errors="coerce")
    meta = (meta.merge(info, on="sample_number", how="left")
            .merge(traits[tcols].drop_duplicates("sample_number"), on="sample_number", how="left")
            .merge(elem[ecols].drop_duplicates("sample_number"), on="sample_number", how="left")
            .merge(fiber[fcols].drop_duplicates("sample_number"), on="sample_number", how="left"))
    # ~220 spectra samples are missing from the field sample list, but the chemistry and
    # trait files carry species_or_type for most of them (chemistry first: where both are
    # in the sample list, the chemistry files agree with it and the trait file sometimes not)
    alt = (pd.concat([d[["sample_number", "species_or_type"]] for d in (elem, fiber, traits)])
           .dropna().drop_duplicates("sample_number").set_index("sample_number")["species_or_type"])
    meta["species_or_type"] = meta["species_or_type"].fillna(meta["sample_number"].map(alt))
    sp = meta["species_or_type"].fillna("").str.split(" - ").str[0].str.strip()
    meta["species_name"] = sp.where(sp.str.match(r"^[A-Z][a-z]+ [a-z]"), None)
    meta["site"] = meta["plot_name"].fillna("").str[:3].map(
        {"DPA": "Dangermond", "DPM": "Dangermond", "DP": "Dangermond", "SRC": "Sedgwick",
         "SRA": "Sedgwick", "CSM": "Carpinteria"}).fillna("Santa Barbara County")
    dinfo = {"title": "SHIFT: Reflectance Measurements for Dried and Ground Leaf Materials; "
                      "SHIFT: Photosynthetic and Leaf Traits; SHIFT: Laboratory Foliar Chemical Analysis",
             "url": "https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=2244",
             "doi": "10.3334/ORNLDAAC/2244; 10.3334/ORNLDAAC/2233; 10.3334/ORNLDAAC/2337",
             "citation": cfg.get("citation"), "license": "NASA open data"}
    return to_tables(cfg, meta, wl, vals, dinfo)
