# Dataset adapter configs

There is one YAML file per dataset; its name is the `dataset_id`. A source loader
(`leafrtm.db.adapters.<source>`) reads the raw files into one metadata row
per spectrum. `core.to_tables` then applies the config. Anything not listed
below is ignored.

```yaml
dataset_id: lopex                 # required; equals the file name
source: ecosis                    # ecosis | spectra_db | ecostress | shift
source_id: <package name or list of packages>
include: true                     # false = documented exclusion (give `exclude_reason`)
title / citation / doi / url / license: optional overrides of source metadata
columns: {plant: "ceil(`Sample_#` / 5)"}  # derived metadata columns (pandas eval, applied first)
filter: ["`Target Type` == 'leaf'"]   # pandas query strings on metadata
sentinels: [-999, -9999]          # values treated as missing (default shown + -99999)

spectra:
  quantity: R | {column: Measurement_type, map: {Reflectance: R, Transmittance: T}}
  scale: auto | percent | fraction   # auto: per spectrum, median > 1.5 => percent
  foreoptic: integrating_sphere | leaf_clip | contact_probe | bare_fiber | other | unknown
  instrument: asd-fs3                # id in vocab/instruments.yaml (or {column, map})
  sample_state: fresh | dried_intact | pressed | ground   (or {column, map})
  target_config: single_leaf | leaf_stack | needle_array | single_needle | unknown
  wavelength_range: [350, 2500]      # optional crop (applied after scaling)
  wavelength_scale: 1000             # multiply wavelengths (e.g. micrometres -> nm)
  drop_integer_bands: true | [0]     # keep only non-integer wavelength columns (all packages,
                                     # or only the listed package indices)
  segments: [1000, 1800]             # detector joins: wavelengths that start a new segment

ids:
  entity: [Sample_ID]     # columns that identify one physical leaf / sample (required)
  time: [Date]            # columns that distinguish repeated measurements of that leaf
  plant: [Tree_ID]        # optional individual-plant id
# Rows sharing entity + time and quantity are replicate spectra of one observation.

species: {name: "Latin Name"} | {name: [Genus, Species]} | {code: "USDA Symbol"} | {constant: "Zea mays"}
         # optional `map: {raw: corrected}` and `code_map: {code: name}`
date: {column: Date, format: "%Y%m%d"} | {year: Year} | 2014 | "2014-07-01"
site: {column: Site} | {name: "Ispra"}
plot: {column: Plot}
lat: Latitude    lon: Longitude    crs: "EPSG:32616"   # crs only for projected coordinates
ecosystem: forest | grassland | ... (schema.ECOSYSTEM)
canopy_position: {column: X, map: {T: top, M: middle, B: bottom}}
sun_shade:       {column: X, map: {...}}
leaf_age:        {column: X, map: {...}}
treatments:
  - {factor: water, column: Treatment, control: ["Well-watered"]}
traits:
  - {column: "LMA (g/m2)", trait: lma, unit: "g/m**2"}       # unit in pint syntax
  - {column: "N_area", trait: n_area, unit: "g/m**2", level: plant}   # level: leaf (default) | plant | plot | species
  - {column: "x", trait: starch_mass, unit: "percent", factor: 0.1}  # optional multiplier before conversion
  - {column: "CHL", trait: chl_area, unit: "mg/m**2", method: ccm300}  # method: non-standard measurement (optical meter etc.)
```

Mapped values (`map`) that aren't listed become missing. Enum values are
validated against `leafrtm.db.schema`. Trait names come from
`vocab/traits.yaml`.

## spectra_db datasets not used

`spectra_db` (pinned commit in `sources/spectra_db.py`) is used only for datasets that EcoSIS
doesn't provide. These directories duplicate EcoSIS packages and are read from EcoSIS
instead:

| spectra_db | EcoSIS dataset(s) used instead |
|---|---|
| `lopex`, `angers` | `lopex`, `angers` |
| `nasa_fft` | `nasa_fft` |
| `nasa_hyspiri` | `uw_bnl_hyspiri_2013_2015` |
| `ngee_arctic` | `ngee_arctic_*` (EcoSIS adds transmittance; only chlorophyll/carotenoids of 35 leaves from 2015 are not on EcoSIS) |
| `ngee_tropics` | `ngee_tropics_panama_2016` (EcoSIS adds PSR+ spectra, LMA, water, C/N) |
| `ecosis_californiatraits` | `hyspiri_california_2013` |
| `ecosis_cedarcreek_biodiversity` | `cedar_creek_esr_2014_contact` |
| `ecosis_cornvarieties` | `corn_varieties_wmars_2014` |
| `ecosis_milkweed_stress` | `milkweed_water_temperature` |
| `ecosis_missoulats` | `missoula_lodgepole_sagebrush_2012` |
| `ecosis_pepper` | `herrmann_pepper_nk` |
| `ecosis_pvy_solanum` | `pvy_potato_hawaii_2014` |
| `ecosis_soybean_aphid` | `soybean_aphid_2013` |
| `ecosis_spectral_variation_leafcanopy` | `amaranth_leaf_canopy_2013` |
| `ecosis_hawaii2000` | excluded (`hawaii_2000_vegetation`: branch/canopy-level spectra) |
| `ecosis_santamonica` | not used: mixed vegetation/soil/NPV library without traits or leaf-level metadata |
| `foster_beetle` | excluded (continuum-removed spectra only) |
