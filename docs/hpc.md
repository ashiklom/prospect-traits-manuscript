# Running the full fit on Yale Bouchet (SLURM + dSQ)

The pilot runs end to end on a workstation (`snakemake --cores 4`). The full fit, every
observation × RTM, is thousands of core-hours, so it runs on Bouchet. The fitting jobs
are submitted with YCRC's [Dead Simple Queue](https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/dsq/).
Everything else runs through Snakemake.

## 1. Set up (login node)

```bash
git clone https://github.com/ashiklom/prospect-traits-manuscript && cd prospect-traits-manuscript
module load miniconda   # or any recent Python; uv manages the rest
curl -LsSf https://astral.sh/uv/install.sh | sh
uv sync                 # creates .venv from uv.lock
```

Downloads need network access, so run them on the login or transfer node. NASA SHIFT needs
`EARTHDATA_LOGIN_USER` and `EARTHDATA_LOGIN_PASSWORD`, or a `~/.netrc` entry for
`urs.earthdata.nasa.gov`.

```bash
.venv/bin/snakemake --cores 2 --config fitset=full -- results/full/fit_set.parquet
```

This runs download → harmonize → selection and writes `results/full/inputs/<batch>.pkl`.

## 2. The full run's configuration

`config/config.yaml` → `fitsets.full` is already configured from the pilot (the comments
there give the evidence):

- **Error model:** heteroskedastic + Matérn-3/2 (`ss`) for every RTM, with correlated R/T
  errors where a leaf has both.
- **Models:** PROSPECT-4, 5, 5B, D and PRO, LIBERTY (Eq. 16 corrected) and LeafSIP.
- **Quantities:** reflectance, transmittance, absorptance and pseudo-absorbance.
- **R only:** a PROSPECT-D fit to reflectance alone on every leaf with R and T.
- **No stage 1:** all runs are stage 2. The Matérn length scale and nugget come from the
  pilot (`discrepancy_from: pilot`). PROSPECT versions not in the pilot borrow
  PROSPECT-D's (`discrepancy.borrow`).

That needs the pilot's learned error structure, `results/pilot/discrepancy/` (7 MB, not in
git). Copy it from the machine that ran the pilot, or unpack the release asset:

```bash
mkdir -p results/pilot && tar -xzf pilot-discrepancy.tar.gz -C results/pilot
```

`n_batches: 400` puts about 230 observations in each job.

## 3. Fits (dSQ)

```bash
.venv/bin/python -m leafrtm.hpc --fitset full --stage 2 --pilot-dir results/pilot
module load dSQ
dsq --job-file hpc/fit_jobs_full_stage2.txt --mem-per-cpu 4g -t 4:00:00 --max-jobs 400 \
    -o hpc/logs/%A_%a.out --batch-file hpc/fits.sh
sbatch hpc/fits.sh
```

`leafrtm.hpc` writes one job line per (run, batch) whose outputs are missing, and prints the
projected cost from the pilot's timings. Each line runs on one core. In October 2026 the
projection was 3,200 jobs and about 3,100 core-hours, roughly 1 core-hour per job. Its
inputs were:

- 91,501 observations, 4,434 of them with R and T;
- the pilot's per-fit times, scaled to the full run's 2 chains × 1,000 iterations.

Per-observation results are cached in `results/full/cache/<run>/`, so a job that hits its
time limit resumes where it stopped.

After the array finishes:

```bash
dsqa -j <jobid> -f hpc/fit_jobs_full_stage2.txt > hpc/rerun.txt   # failed or timed-out jobs
# resubmit hpc/rerun.txt the same way, until it is empty
```

To add stage-1 (OU) runs and learn the error structure from the full data instead, give the
runs `stage: 1`, drop `discrepancy_from`, and run
`snakemake --config fitset=full -- results/full/discrepancy/report.json` between the two job
files (`--stage 1`, then `--stage 2`).

## 4. Analyses and site

The remaining steps (aggregate → analysis → site_data → site) are light:

```bash
.venv/bin/snakemake --cores 4 --config fitset=full
```

Publish the site data as a release asset (`results/full/site-data-full.tar.gz`; see the
README). Then render the site anywhere with `LEAFRTM_FITSET=full`.
