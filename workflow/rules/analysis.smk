# Analyses, site data and the website.

ANALYSIS_TABLES = ["estimates", "validation_leaf", "validation_metrics", "validation_species",
                   "spectral_residuals", "heldout", "data_use_comparison", "intercompare",
                   "diagnostics", "error_model_comparison", "best_config", "best_estimates",
                   "variance", "treatments", "correlations", "within_species_correlations"]


rule analysis:
    input:
        "results/{fitset}/summary/info.parquet",
        ancient("data/processed/db/report.json"),
    output:
        expand("results/{{fitset}}/analysis/{t}.parquet", t=ANALYSIS_TABLES),
    shell:
        "{PY} -m leafrtm.analysis.compute --fitset-dir results/{wildcards.fitset}"


rule recovery:
    """Synthetic recovery of known parameters (pilot only)."""
    output:
        "results/pilot/analysis/recovery.parquet",
    threads: 1
    shell:
        "OMP_NUM_THREADS=1 {PY} -m leafrtm.analysis.recovery --out results/pilot/analysis"


rule feasibility:
    """Absorptance / pseudo-absorbance feasibility (fit set `feasibility`)."""
    input:
        "results/feasibility/summary/info.parquet",
    output:
        "results/pilot/analysis/apa_feasibility.parquet",
    shell:
        "{PY} -m leafrtm.analysis.feasibility --fitset-dir results/feasibility --out results/pilot/analysis"


def site_data_inputs(wc):
    deps = [f"results/{wc.fitset}/analysis/{t}.parquet" for t in ANALYSIS_TABLES]
    if wc.fitset == "pilot":
        deps += ["results/pilot/analysis/recovery.parquet",
                 "results/pilot/analysis/apa_feasibility.parquet"]
    return deps


rule site_data:
    input:
        site_data_inputs,
    output:
        directory("results/{fitset}/site_data"),
        "results/{fitset}/site-data-{fitset}.tar.gz",
    shell:
        "{PY} -m leafrtm.site_data --fitset-dir results/{wildcards.fitset}"


rule site:
    input:
        "results/{fitset}/site_data",
        glob.glob("site/*.qmd"),
        "site/_lib.py",
        "site/_quarto.yml",
    output:
        "results/{fitset}/site/index.html",
    shell:
        "cd site && LEAFRTM_FITSET={wildcards.fitset} LEAFRTM_SITE_DATA=../results/{wildcards.fitset}/site_data "
        "QUARTO_PYTHON=$(which python) quarto render --output-dir ../results/{wildcards.fitset}/site"
