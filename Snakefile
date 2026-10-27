# Leaf RTM trait-retrieval workflow.
#
#   snakemake --cores 4                       # pilot (default fitset), through the website
#   snakemake --cores 4 --config fitset=full  # full run (see docs/hpc.md for Bouchet/dSQ)
#   snakemake -n                              # dry run
#
# data/ (raw, interim, processed) and results/ are not version controlled.

configfile: "config/config.yaml"

import sys
from pathlib import Path

FITSET = config["fitset"]
# Rules run in the interpreter snakemake itself runs in (e.g. `uv run snakemake`), activated or not.
PYEXE = sys.executable
PYBIN = str(Path(PYEXE).parent)
PY = f"{PYEXE} -W ignore"


def fs(wc):
    return config["fitsets"][wc.fitset]


def runs(fitset, stage=None, **match):
    out = []
    for r in config["fitsets"][fitset]["runs"]:
        if stage is not None and r["stage"] != stage:
            continue
        if all(r.get(k) == v for k, v in match.items()):
            out.append(r["id"])
    return out


def batches(fitset):
    return range(config["fitsets"][fitset]["n_batches"])


include: "workflow/rules/data.smk"
include: "workflow/rules/fit.smk"
include: "workflow/rules/analysis.smk"


rule all:
    input:
        f"results/{FITSET}/site/index.html",
