"""Yang et al. (2016) leaf spectra and traits, from the authors' data repository
(github.com/ashiklom/LeafTraitSpectraSeasonality-Yangetal2016RSE), pinned to a commit.

Files are fetched from raw.githubusercontent.com (no credentials) into
``data/raw/yang2016/`` and checked against SHA-256 sums; the two zip archives are
extracted to ``2011-MV/`` and ``2012_HF/`` (macOS ``__MACOSX`` metadata skipped).
"""

from __future__ import annotations

import hashlib
import shutil
import zipfile
from pathlib import Path, PurePosixPath

from leafrtm.db.sources.common import download, raw_dir

REPO = "ashiklom/LeafTraitSpectraSeasonality-Yangetal2016RSE"
COMMIT = "97312e0f4ef4a942d716938793db4c0b01f5f162"
BASE = f"https://raw.githubusercontent.com/{REPO}/{COMMIT}"
FILES = {  # name -> sha256 at COMMIT
    "2011-MV.zip": "c3e51aff9d6ad9a668e9cafd2ff869bffead531698a86ff4f4742f9fae1c8299",
    "2012_HF.zip": "73d17544d716aa3664adea98eb05be4fce513221811d744ea527bd2a3e2b8d21",
    "2011_MV_leaftraits.xlsx": "90bf7374ebcbb2823536d4b944e29debc7d87bcff3344c4a4b32c70ffedbcd85",
    "2012_HF_leaftraits.xlsx": "08c8bc51c4cdcc6bfa1f7bb99934903f7241033bb78795fc448010830d2538d1",
    "2012HFReadMe.txt": None,  # text file (line endings may be normalized); not checked
    "README.md": None,
}
ZIPS = {"2011-MV.zip": "2011-MV", "2012_HF.zip": "2012_HF"}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _extract(zip_path: Path, dest: Path) -> None:
    """Extract to ``dest`` (atomically, via a temporary directory)."""
    tmp = dest.with_name(dest.name + ".part")
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    with zipfile.ZipFile(zip_path) as z:
        for info in z.infolist():
            p = PurePosixPath(info.filename)
            if info.is_dir() or p.parts[0] == "__MACOSX" or p.name.startswith("._"):
                continue
            if p.is_absolute() or ".." in p.parts:
                raise ValueError(f"unsafe path in {zip_path.name}: {info.filename}")
            out = tmp.joinpath(*p.parts)
            out.parent.mkdir(parents=True, exist_ok=True)
            with z.open(info) as src, open(out, "wb") as dst:
                shutil.copyfileobj(src, dst)
    tmp.rename(dest)


def fetch_all() -> Path:
    d = raw_dir("yang2016")
    for name, sha in FILES.items():
        p = download(f"{BASE}/{name}", d / name)
        if sha and _sha256(p) != sha:
            p.unlink()
            raise ValueError(f"checksum mismatch for {name} (expected commit {COMMIT})")
    for name, sub in ZIPS.items():
        if not (d / sub).exists():
            _extract(d / name, d / sub)
    return d
