"""Summarize an EcoSIS export's metadata columns (for writing adapter configs).

Usage: python scripts/inspect_ecosis.py <package-name> [max_values]
"""

import sys

import pandas as pd

from leafrtm.db.sources.common import raw_dir


def is_wavelength(c):
    try:
        float(c)
        return True
    except ValueError:
        return False


def main(name, k=8):
    df = pd.read_csv(raw_dir("ecosis") / f"{name}.csv.gz", low_memory=False)
    wl = [c for c in df.columns if is_wavelength(c)]
    meta = [c for c in df.columns if not is_wavelength(c)]
    w = pd.Series([float(c) for c in wl])
    print(f"# {name}\nrows={len(df)} bands={len(wl)} range={w.min()}-{w.max()} "
          f"step(median)={w.diff().median()} duplicates={w.duplicated().sum()}")
    v = df[wl].apply(pd.to_numeric, errors="coerce")
    print(f"value range: {v.min().min():.4g} .. {v.max().max():.4g}; median {v.stack().median():.4g}")
    for c in meta:
        s = df[c]
        u = s.dropna().unique()
        ex = ", ".join(map(str, u[:k]))
        print(f"- {c!r}: n_unique={len(u)} n_na={s.isna().sum()} e.g. {ex[:200]}")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 8)
