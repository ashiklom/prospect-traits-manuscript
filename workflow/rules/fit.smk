# Inversions. Each job fits one batch for one run (hetero + OU errors; see
# leafrtm.fit); per-observation results are cached under
# results/<fitset>/cache/<run>/ so interrupted jobs resume.


def fit_inputs(wc):
    return [ancient(f"results/{wc.fitset}/fit_set.parquet")]


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


rule aggregate:
    input:
        lambda wc: fits_of(wc.fitset, runs(wc.fitset)),
    output:
        expand("results/{{fitset}}/summary/{part}.parquet", part=["summary", "draws", "spectra", "info"]),
    shell:
        "{PY} -m leafrtm.aggregate --fitset-dir results/{wildcards.fitset}"
