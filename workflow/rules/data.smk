# Download, harmonize and select.
import glob

ADAPTER_CONFIGS = sorted(glob.glob("src/leafrtm/db/adapters/configs/*.yaml"))


rule harmonize:
    """Download every source (cached), convert each dataset and merge."""
    input:
        ADAPTER_CONFIGS,
        "src/leafrtm/db/vocab/traits.yaml",
    output:
        "data/processed/db/report.json",
    threads: 4
    shell:
        "{PY} -m leafrtm.db.harmonize"


rule select:
    input:
        ancient("data/processed/db/report.json"),
    output:
        fit_set="results/{fitset}/fit_set.parquet",
    shell:
        "{PY} -m leafrtm.select --fitset {wildcards.fitset}"
