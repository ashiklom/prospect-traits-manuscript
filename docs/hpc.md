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

## 2. Configure the full run

Edit `config/config.yaml` → `fitsets.full.runs`, using the pilot's
[Error models](../site/error-models.qmd) comparison to pick the error model(s), and the R vs
R+T comparison to decide whether R-only fits are needed. For example, for every RTM:

```yaml
- {id: prospectd_ou, model: prospectd, errmodel: ou, rt: false, data_use: all, stage: 1}
- {id: prospectd_ld_rt, model: prospectd, errmodel: ld, rt: true, data_use: all, stage: 2}
```

Stage-1 `ou` runs are needed by `learn` whenever any stage-2 `ss` or `ld` run uses that RTM.
`n_batches` sets the job size: 400 batches puts about 100–150 observations in each.

## 3. Stage 1 fits (dSQ)

```bash
.venv/bin/python -m leafrtm.hpc --fitset full --stage 1 --pilot-dir results/pilot   # prints projected core-hours
module load dSQ
dsq --job-file hpc/fit_jobs_full_stage1.txt --mem-per-cpu 4g -t 12:00:00 --max-jobs 400 \
    -o hpc/logs/%A_%a.out --batch-file hpc/stage1.sh
sbatch hpc/stage1.sh
```

Each job line runs one (run, batch) on one core. Per-observation results are cached in
`results/full/cache/<run>/`, so a job that hits its time limit resumes where it stopped.

After the array finishes:

```bash
dsqa -j <jobid> -f hpc/fit_jobs_full_stage1.txt > hpc/rerun_stage1.txt   # failed or timed-out jobs
# resubmit hpc/rerun_stage1.txt the same way, until it is empty
```

## 4. Learn error-model structure, then stage 2

```bash
.venv/bin/snakemake --cores 4 --config fitset=full -- results/full/discrepancy/report.json
.venv/bin/python -m leafrtm.hpc --fitset full --stage 2
dsq --job-file hpc/fit_jobs_full_stage2.txt --mem-per-cpu 4g -t 12:00:00 --max-jobs 400 \
    -o hpc/logs/%A_%a.out --batch-file hpc/stage2.sh && sbatch hpc/stage2.sh
```

## 5. Analyses and site

The remaining steps (aggregate → analysis → site_data → site) are light:

```bash
.venv/bin/snakemake --cores 4 --config fitset=full
```

Publish the site data as a release asset (`results/full/site-data-full.tar.gz`; see the
README). Then render the site anywhere with `LEAFRTM_FITSET=full`.
