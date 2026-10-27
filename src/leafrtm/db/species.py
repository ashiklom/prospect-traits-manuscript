"""Species name resolution.

Pass 1 (during harmonization): raw names and USDA PLANTS codes are cleaned
to a binomial and matched against the USDA PLANTS checklist shipped with
spectra_db (``species_info.csvy``). Names not found there are listed in
``vocab/species_resolution.csv`` and resolved by hand (manual review).

Pass 2 (``python -m leafrtm.db.species gbif``): every resolved name is
checked against the GBIF backbone (accepted name, family, match type); the
results are committed as ``vocab/species_gbif.csv`` so rebuilds are offline.

Species attributes (growth form, leaf type, phenology, photosynthetic
pathway) come from ``vocab/species_attributes.csv`` (by species, falling back
to genus).
"""

from __future__ import annotations

import re
import sys
import time
from functools import cache
from importlib.resources import files

import pandas as pd
import requests

from leafrtm.db.sources import spectra_db

VOCAB = files("leafrtm.db") / "vocab"
RESOLUTION = "species_resolution.csv"
GBIF = "species_gbif.csv"
ATTRIBUTES = "species_attributes.csv"

_AUTHOR = re.compile(r"\s+(\(?[A-Z][\w.\-'& ]*\)?\s*)+$")
_QUOTED = re.compile(r"['\"‘’“”].*?['\"‘’“”]")


def clean_name(raw) -> str | None:
    """Reduce a raw taxon string to 'Genus species' (or 'Genus' for 'Genus sp.')."""
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return None
    s = str(raw).replace("_", " ").replace("×", "x ").strip()
    s = _QUOTED.sub("", s)  # cultivars
    s = re.sub(r"\s+", " ", s)
    toks = s.split(" ")
    if not toks or not re.match(r"^[A-Z][a-z]+$", toks[0].capitalize()):
        return None
    genus = toks[0].capitalize()
    if len(toks) == 1 or toks[1].lower() in ("sp", "sp.", "spp", "spp.", "species", "cf.", "aff."):
        return genus
    ep = toks[1].lower()
    if ep == "x" and len(toks) > 2:  # hybrid: Genus x epithet
        ep = "x " + toks[2].lower()
    if not re.match(r"^(x )?[a-z][a-z\-]+$", ep):
        return genus
    return f"{genus} {ep}"


def species_id(name: str | None) -> str | None:
    if not isinstance(name, str) or not name:
        return None
    return name.lower().replace(" ", "_")


@cache
def usda() -> pd.DataFrame:
    """USDA PLANTS checklist from spectra_db (code, scientific name, family)."""
    from leafrtm.db.adapters.spectra_db import read_csvy

    _, df = read_csvy(spectra_db.checkout() / "species_info.csvy")
    df = df[["species_code", "scientific_name", "family", "Database_ID"]].copy()
    df["name"] = df["scientific_name"].map(clean_name)
    return df.dropna(subset=["name"])


def _read_vocab(name, **kw) -> pd.DataFrame:
    p = VOCAB / name
    if not p.is_file():
        return pd.DataFrame()
    with p.open() as f:
        return pd.read_csv(f, **kw)


def resolve(raw: pd.DataFrame) -> pd.DataFrame:
    """Map raw (name, code) pairs to resolved names.

    Returns raw with columns ``resolved`` and ``resolution`` (usda_code,
    usda_name, manual, cleaned, unresolved)."""
    u = usda()
    by_code = dict(zip(u["species_code"].str.upper(), u["name"]))
    usda_names = set(u["name"])
    manual = _read_vocab(RESOLUTION, dtype=str, keep_default_na=False)
    man = {}
    for r in manual.itertuples():
        man[(r.raw_name, r.raw_code)] = r.resolved_name or None
    out = raw.copy()
    if "name" not in out:
        out["name"] = None
    if "code" not in out:
        out["code"] = None
    res, how = [], []
    for name, code in zip(out["name"], out["code"]):
        key = ("" if pd.isna(name) or name is None else str(name),
               "" if pd.isna(code) or code is None else str(code))
        if key in man:
            res.append(man[key])
            how.append("manual")
            continue
        if key[1] and key[1].upper() in by_code:
            res.append(by_code[key[1].upper()])
            how.append("usda_code")
            continue
        c = clean_name(key[0])
        if c in usda_names:
            res.append(c)
            how.append("usda_name")
        elif c:
            res.append(c)
            how.append("cleaned")  # not in USDA: needs manual review
        else:
            res.append(None)
            how.append("unresolved")
    out["resolved"] = res
    out["resolution"] = how
    return out


def review_table(resolved: pd.DataFrame) -> pd.DataFrame:
    """Rows that need manual review (not matched to USDA PLANTS)."""
    r = resolved[~resolved["resolution"].isin(["usda_code", "usda_name", "manual"])]
    return (r.groupby(["name", "code"], dropna=False)
            .agg(datasets=("dataset_id", lambda s: ";".join(sorted(set(s)))),
                 n_entities=("entity_id", "size"), suggestion=("resolved", "first"))
            .reset_index().sort_values("n_entities", ascending=False))


# ---------------------------------------------------------------------------
# Pass 2: GBIF backbone
# ---------------------------------------------------------------------------


def gbif_match(name: str, session=None) -> dict:
    s = session or requests
    r = s.get("https://api.gbif.org/v1/species/match",
              params={"name": name, "kingdom": "Plantae", "strict": "false"}, timeout=60)
    r.raise_for_status()
    d = r.json()
    accepted = d.get("species") if d.get("rank") in ("SPECIES", "SUBSPECIES", "VARIETY", "FORM") \
        else d.get("canonicalName")
    return {"name": name, "gbif_accepted": accepted, "gbif_status": d.get("status"),
            "gbif_match": d.get("matchType"), "gbif_confidence": d.get("confidence"),
            "gbif_rank": d.get("rank"), "family": d.get("family"), "genus": d.get("genus"),
            "order": d.get("order"), "gbif_key": d.get("usageKey")}


def gbif_table(names, existing: pd.DataFrame | None = None) -> pd.DataFrame:
    done = set(existing["name"]) if existing is not None and len(existing) else set()
    s = requests.Session()
    rows = []
    for n in sorted(set(names) - done):
        rows.append(gbif_match(n, s))
        time.sleep(0.05)
    new = pd.DataFrame(rows)
    return pd.concat([existing, new], ignore_index=True) if existing is not None else new


def main(argv=None):
    argv = argv or sys.argv[1:]
    if argv and argv[0] == "gbif":
        from leafrtm.db.harmonize import processed_dir
        sp = pd.read_parquet(processed_dir() / "species.parquet")
        old = _read_vocab(GBIF)
        tab = gbif_table(sp["scientific_name"], old if len(old) else None)
        tab.sort_values("name").to_csv(sys.stdout if len(argv) < 2 else argv[1], index=False)


if __name__ == "__main__":
    main()
