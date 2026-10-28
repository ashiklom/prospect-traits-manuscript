# Inversions. Each job fits one batch for one run; per-observation results are
# cached under results/<fitset>/cache/<run>/ so interrupted jobs resume.


def fit_inputs(wc):
    deps = [ancient(f"results/{wc.fitset}/fit_set.parquet")]
    stage = {r["id"]: r["stage"] for r in config["fitsets"][wc.fitset]["runs"]}[wc.run]
    if stage == 2:  # learned error structure, possibly from another fit set (discrepancy_from)
        src = config["fitsets"][wc.fitset].get("discrepancy_from", wc.fitset)
        deps.append(f"results/{src}/discrepancy/report.json" if src == wc.fitset
                    else ancient(f"results/{src}/discrepancy/report.json"))
    return deps


rule fit:
    input:
        fit_inputs,
    output:
        multiext("results/{fitset}/fits/{run}/batch{batch}", ".summary.parquet",
                 ".draws.parquet", ".spectra.parquet", ".info.parquet"),
    wildcard_constraints:
        batch=r"\d{4}",
    threads: 1
    resources:
        mem_mb=3000,
    shell:
        "OMP_NUM_THREADS=1 {PY} -m leafrtm.run --fitset-dir results/{wildcards.fitset} "
        "--run {wildcards.run} --batch {wildcards.batch}"


def fits_of(fitset, run_ids):
    return expand("results/{fitset}/fits/{run}/batch{batch:04d}.summary.parquet",
                  fitset=fitset, run=run_ids, batch=batches(fitset))


def stage1_learn_runs(fitset):
    return runs(fitset, stage=1, errmodel="ou", rt=False, data_use="all")


rule learn:
    """Error-model structure (SS correlation, LD discrepancies) from stage-1 ou fits."""
    input:
        lambda wc: fits_of(wc.fitset, stage1_learn_runs(wc.fitset)),
    output:
        "results/{fitset}/discrepancy/report.json",
        "results/{fitset}/discrepancy/ss_correlation.json",
    params:
        runs=lambda wc: " ".join(stage1_learn_runs(wc.fitset)),
    shell:
        "{PY} -m leafrtm.learn --fitset-dir results/{wildcards.fitset} --stage1-runs {params.runs}"


rule aggregate:
    input:
        lambda wc: fits_of(wc.fitset, runs(wc.fitset)),
    output:
        expand("results/{{fitset}}/summary/{part}.parquet", part=["summary", "draws", "spectra", "info"]),
    shell:
        "{PY} -m leafrtm.aggregate --fitset-dir results/{wildcards.fitset}"
