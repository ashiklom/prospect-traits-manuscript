"""Assemble the tables the website needs into one directory (+ a tarball).

    python -m leafrtm.site_data --fitset-dir results/pilot

Copies ``<fitset>/analysis/*.parquet``, adds database inventory tables and a
sample of example spectral fits, and writes ``site-data-<fitset>.tar.gz`` for
publishing as a GitHub release asset.
"""

from __future__ import annotations

import argparse
import shutil
import tarfile
from pathlib import Path

import pandas as pd

from leafrtm.analysis.inventory import main as inventory_main

N_EXAMPLE_LEAVES = 40


def discrepancy_table(root: Path) -> pd.DataFrame:
    """Learned discrepancies (first available fold of each model/tag/class) in long form."""
    from leafrtm.discrepancy import Discrepancy

    rows = []
    first = {}
    for p in sorted(root.glob("*/*/*/fold*.npz")):
        first.setdefault(p.parent, p)
    for p in first.values():
        model, tag, klass = p.parts[-4:-1]
        dd = Discrepancy.load(p)
        ng = len(dd.grid)
        for qi, q in enumerate(dd.quantities):
            sl = slice(qi * ng, (qi + 1) * ng)
            df = pd.DataFrame({"model": model, "tag": tag, "klass": klass, "key": f"{model}/{tag}/{klass}",
                               "quantity": q, "wavelength": dd.grid, "mu": dd.mu[sl],
                               "sd": dd.D[sl] ** 0.5, "ell": dd.ell, "nugget": dd.nugget,
                               "n_leaves": dd.n_leaves})
            for j in range(dd.Lambda.shape[1]):
                df[f"factor{j + 1}"] = dd.Lambda[sl, j]
            rows.append(df)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def liberty_table(d: Path) -> pd.DataFrame:
    """Per LIBERTY fit: posterior-median parameters and mean observed vs fitted
    reflectance over 800-900 nm (the bands nearest LIBERTY's validity edge)."""
    summ, spec = d / "summary" / "summary.parquet", d / "summary" / "spectra.parquet"
    if not (summ.exists() and spec.exists()):
        return pd.DataFrame()
    sm = pd.read_parquet(summ, filters=[("model", "==", "liberty"), ("is_rtm", "==", True)])
    if sm.empty:
        return pd.DataFrame()
    par = sm.pivot_table(index=["observation_id", "run_id"], columns="param", values="q50").reset_index()
    sp = pd.read_parquet(spec, filters=[("model", "==", "liberty"), ("quantity", "==", "R")])
    sp = sp[sp["wavelength"].between(800, 900)]
    nir = (sp.groupby(["observation_id", "run_id"])[["y", "yhat_q50"]].mean()
           .rename(columns={"y": "nir_obs", "yhat_q50": "nir_fit"}).reset_index())
    return par.merge(nir, on=["observation_id", "run_id"], how="left")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fitset-dir", required=True)
    a = ap.parse_args(argv)
    d = Path(a.fitset_dir)
    out = d / "site_data"
    out.mkdir(parents=True, exist_ok=True)
    for p in (d / "analysis").glob("*.parquet"):
        shutil.copy(p, out / p.name)
    inventory_main(["--out", str(out)])
    spec = d / "summary" / "spectra.parquet"
    if spec.exists():
        s = pd.read_parquet(spec)
        if len(s):
            leaves = (s["observation_id"].drop_duplicates()
                      .sample(min(N_EXAMPLE_LEAVES, s["observation_id"].nunique()), random_state=1))
            s[s["observation_id"].isin(leaves)].to_parquet(out / "spectra_examples.parquet", index=False)
    lib = liberty_table(d)
    if len(lib):
        lib.to_parquet(out / "liberty_fits.parquet", index=False)
    learned = discrepancy_table(d / "discrepancy")
    if len(learned):
        learned.to_parquet(out / "discrepancy_learned.parquet", index=False)
    tar = d / f"site-data-{d.name}.tar.gz"
    with tarfile.open(tar, "w:gz") as t:
        t.add(out, arcname="site_data")
    print(f"site data -> {out} and {tar}")


if __name__ == "__main__":
    main()
