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

- **Errors:** heteroskedastic + Ornstein–Uhlenbeck, independently for R, T, A and PA. This
  is the only error model fitted after commit `0cc1e51`.
- **Models:** PROSPECT-4, 5, 5B, D and PRO, LIBERTY (Eq. 16 corrected) and LeafSIP, fitted
  to every observation.
- **R only:** a PROSPECT-D fit to reflectance alone on every leaf with R and T.

All runs are independent, so there is a single job file. `n_batches: 400` puts
about 230 observations in each job.

## 3. Fits (dSQ)

```bash
.venv/bin/python -m leafrtm.hpc --fitset full --pilot-dir results/pilot
module load dSQ
dsq --job-file hpc/fit_jobs_full.txt --mem-per-cpu 4g -t 4:00:00 --max-jobs 400 \
    -o hpc/logs/%A_%a.out --batch-file hpc/fits.sh
sbatch hpc/fits.sh
```

`leafrtm.hpc` writes one job line per (run, batch) whose outputs are missing. It also
prints the projected cost from the pilot's OU timings, scaled to the full run's 2 chains ×
1,000 iterations. The projection needs `results/pilot/summary/info.parquet`; drop
`--pilot-dir` without it. In October 2026 it was:

- 91,501 observations, 4,434 of them with R and T;
- 3,200 jobs and about 3,400 core-hours, roughly 1 core-hour per job.

Per-observation results are cached in `results/full/cache/<run>/`, so a job that hits its
time limit resumes where it stopped.

After the array finishes:

```bash
dsqa -j <jobid> -f hpc/fit_jobs_full.txt > hpc/rerun.txt   # failed or timed-out jobs
# resubmit hpc/rerun.txt the same way, until it is empty
```

### Earlier error models

Up to commit `0cc1e51` (tag `v0.3-pilot-all-errmodels`), the fitting code could also fit
three alternatives:

- a Matérn-3/2 state-space error model;
- a learned shared discrepancy;
- correlated R/T errors.

That needed two stages: OU fits, then `leafrtm.learn`, then the alternatives. Their
likelihoods and the learning code remain in `likelihood.py`, `discrepancy.py` and
`learn.py` as a reference. To run them, check out that commit.

Only the branch was pushed from the session that made the tag. If the tag is missing on
GitHub, recreate it:

```bash
git tag -a v0.3-pilot-all-errmodels 0cc1e51 -m "Last version that dispatches on error model"
git push origin v0.3-pilot-all-errmodels
```

## 4. Analyses and site

The remaining steps (aggregate → analysis → site_data → site) are light:

```bash
.venv/bin/snakemake --cores 4 --config fitset=full
```

Publish the site data as a release asset (`results/full/site-data-full.tar.gz`; see the
README). Then render the site anywhere with `LEAFRTM_FITSET=full`.
