# Chloroplasts to canopies: Leaf optical properties shed light on foliar trait variability at individual to global scales

Source code and ancillary files for the manuscript
"Chloroplasts to canopies: Leaf optical properties shed light on foliar trait variability at
individual to global scales" by Alexey N. Shiklomanov and colleagues.

The repository contains three things:

1. **A harmonized database** of leaf spectra and paired trait measurements, built entirely by
   code from EcoSIS, the ECOSTRESS spectral library, NASA SHIFT (ORNL DAAC) and
   [`spectra_db`](https://github.com/ashiklom/spectra_db).
2. **Bayesian inversions** of seven leaf radiative transfer models (RTMs): PROSPECT-4, 5,
   5B, D and PRO, LIBERTY and LeafSIP. They use PyMC with the nutpie sampler, every replicate
   R and T spectrum of a leaf, and structured error models.
3. **A static website** with interactive analyses, plus the manuscript (Quarto).

Everything from download to website is orchestrated with [Snakemake](https://snakemake.readthedocs.io).
Raw data, processed data and fit results are **not** version controlled. The tables behind
the website are published as a GitHub release asset.

## Quick start

```bash
uv sync                                   # Python 3.13 environment from uv.lock
uv run pytest                             # forward models, likelihoods, inversion, observation operator
uv run snakemake --cores 4                # pilot: download -> harmonize -> fit -> analyses -> website
uv run snakemake --cores 4 --config fitset=smoke   # tiny end-to-end test
```

* **Website output:** `results/<fitset>/site/index.html`.
* **Credentials:** NASA SHIFT downloads need Earthdata login credentials in
  `EARTHDATA_LOGIN_USER`/`EARTHDATA_LOGIN_PASSWORD` (or `~/.netrc`).
* **Full run:** it is HPC-sized; see [`docs/hpc.md`](docs/hpc.md) for the Yale Bouchet / dSQ
  workflow.

To render the website from published results instead of re-running the fits:

```bash
gh release download <tag> -p 'site-data-pilot.tar.gz'
mkdir -p results/pilot && tar -xzf site-data-pilot.tar.gz -C results/pilot
uv run snakemake --cores 1 --allowed-rules site -- results/pilot/site/index.html
```

## Why Snakemake?

* **It is Python-native.** Rules call the package directly; there is no server or scheduler
  service, unlike Airflow, Prefect or Dagster.
* **Its file-based DAG fits the work.** Wildcards express dataset → harmonize →
  (RTM run × batch) → summarize → website.
* **It resumes cleanly.** Interrupted runs pick up where they stopped, which matters in
  ephemeral containers and on clusters.
* **It scales.** The same rules run on a laptop and on a SLURM cluster; the per-batch fit
  commands are also emitted as dSQ job files.

`workflow/profiles/default` restricts rerun triggers to inputs and timestamps, so editing
code doesn't silently trigger thousands of refits.

## Layout

```
Snakefile, workflow/rules/          workflow
config/config.yaml                  fit sets (smoke, pilot, feasibility, full), runs, sampler settings
config/priors.yaml                  priors for RTM and error-model parameters
config/rtm_traits.yaml              RTM parameter <-> measured trait mapping (units)
src/leafrtm/rtm/                    PROSPECT 4/5/5B/D/PRO, LIBERTY, LeafSIP (NumPy + PyTensor); data/PROVENANCE.md
src/leafrtm/obsop.py                observation operator (model grid -> observed bands; no resampling of data)
src/leafrtm/likelihood.py           error models: OU, Matern state-space (Vecchia), learned discrepancy (Woodbury)
src/leafrtm/discrepancy.py          learning shared RTM discrepancy from first-stage residuals
src/leafrtm/fit.py, run.py          PyMC model builder, nutpie inverter (compile once, swap data per leaf)
src/leafrtm/db/                     sources (downloaders), adapters (+ one YAML config per dataset), schema, harmonize, QC, species
src/leafrtm/analysis/               analyses for the website (objectives A and B)
site/                               Quarto website
docs/hpc.md                         running the full fit on Bouchet
latex/, prospect-traits-manuscript.qmd   manuscript
```

## Database organization

The schema is defined in `src/leafrtm/db/schema.py`; the tables are Parquet files under
`data/processed/db/`.

* **Hierarchy:**
  * `dataset` → `site` → `plot`
  * `entity` (a physical leaf or sample: species, plot, plant) → `observation` (an entity
    at one time and sample state; treatments)
  * `observation` → `spectra` (replicate measurements, each on its own native wavelength
    grid) and `traits` (canonical units; leaf or pooled-sample level)
* **IDs:** `spectra_id` → `observation_id` → `entity_id`, formatted `<dataset_id>:<local id>`.
  RTMs are inverted per `observation_id`, using all of its spectra.
* **Enums:** categorical metadata (sample state, target configuration, foreoptic, quantity,
  canopy position, treatment factor, growth form, …) use fixed enumerations.

## Manuscript

Run `quarto render` in the repository root to build PDF, HTML and DOCX versions of the
manuscript in `_output/`.
