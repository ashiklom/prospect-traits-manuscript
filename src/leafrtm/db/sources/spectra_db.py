"""ashiklom/spectra_db: curated leaf spectra datasets in CSVY format.

Only datasets that are not available from EcoSIS are used from here (see the
adapter configs). The repository is cloned at a pinned commit.
"""

from __future__ import annotations

import subprocess

from leafrtm.db.sources.common import raw_dir

REPO = "https://github.com/ashiklom/spectra_db"
COMMIT = "79f8c0213702a79790208829043a49eefe807356"


def checkout():
    d = raw_dir("spectra_db") / "repo"
    if not (d / ".git").exists():
        subprocess.run(["git", "clone", "-q", REPO, str(d)], check=True)
    head = subprocess.run(["git", "-C", str(d), "rev-parse", "HEAD"], check=True,
                          capture_output=True, text=True).stdout.strip()
    if head != COMMIT:
        subprocess.run(["git", "-C", str(d), "fetch", "-q", "origin"], check=True)
        subprocess.run(["git", "-C", str(d), "checkout", "-q", COMMIT], check=True)
    return d
