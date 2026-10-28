"""ECOSTRESS spectral library (speclib.jpl.nasa.gov; Meerdink et al. 2019).

There is no API: the list of spectra comes from the site's search endpoint,
and individual spectrum files from ``/ecospeclibdata/<name>``. Only VSWIR
spectra of the vegetation and non-photosynthetic vegetation classes are
fetched.
"""

from __future__ import annotations

import re
import time

import requests

from leafrtm.db.sources.common import get, raw_dir

BASE = "https://speclib.jpl.nasa.gov"
CLASSES = ("vegetation", "nonphotosyntheticvegetation")


def list_spectra(searchtype: str, session=None) -> list[str]:
    """All VSWIR spectrum file names of one class (the site pages by 100:
    ``maxhits`` is the end of the requested page)."""
    s = session or requests.Session()
    names: set[str] = set()
    page = 100
    while True:
        r = s.post(f"{BASE}/ecospeclib_list", timeout=120, data={
            "searchtype": searchtype, "classsel": "All", "subclass": "All", "mname": "",
            "xstart": "", "xstop": "", "maxhits": str(page), "wavelength": "Any"})
        r.raise_for_status()
        new = set(re.findall(r"getplot\('([^']+)'\)", r.text)) - names
        if not new:
            break
        names |= new
        page += 100
    return sorted(n for n in names if ".vswir." in n)


def fetch_all() -> list[str]:
    s = requests.Session()
    out = []
    for cls in CLASSES:
        d = raw_dir("ecostress") / cls
        d.mkdir(parents=True, exist_ok=True)
        for name in list_spectra(cls, s):
            p = d / name
            if not p.exists():
                p.write_text(get(f"{BASE}/ecospeclibdata/{name}", session=s).text)
                time.sleep(0.2)
            out.append(str(p))
    return out
