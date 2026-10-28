"""NASA SHIFT campaign leaf data (ORNL DAAC; Earthdata login required).

Uses EARTHDATA_LOGIN_USER / EARTHDATA_LOGIN_PASSWORD (or ~/.netrc) for
urs.earthdata.nasa.gov. Datasets: 2244 (dried & ground leaf reflectance),
2233 (LMA, LWC, chlorophyll), 2337 (foliar chemistry), 2295 (plots, samples).
"""

from __future__ import annotations

import os

import requests

from leafrtm.db.sources.common import download, raw_dir

URS = "urs.earthdata.nasa.gov"
BASE = "https://data.ornldaac.earthdata.nasa.gov/protected/shift"
FILES = {
    "SHIFT_DriedGround_Leaf_Reflec": ["Ovendried_Spectra.csv", "Flashfrozen_Spectra.csv"],
    "SHIFT_Leaf_Traits_Chl_SB_CA": ["SHIFT_Leaf_Traits_LMA_LWC_Chl.csv"],
    "SHIFT_Foliar_Chemical_Analysis": ["OD_elemental_daac.csv", "FF_elemental_daac.csv",
                                       "OD_AnkomFiber_daac.csv", "FF_AnkomFiber_daac.csv"],
    "SHIFT_Field_Survey_Metadata": ["SHIFT_vegetation_sample_list.csv",
                                    "SHIFT_vegetation_shift_plots.csv",
                                    "SHIFT_vegetation_species_list.csv",
                                    "SHIFT_vegetation_plot_event_list.csv"],
}


class EarthdataSession(requests.Session):
    """Sends credentials only to the Earthdata login host across redirects."""

    def __init__(self):
        super().__init__()
        user, pw = os.environ.get("EARTHDATA_LOGIN_USER"), os.environ.get("EARTHDATA_LOGIN_PASSWORD")
        if user and pw:
            self.auth = (user, pw)

    def rebuild_auth(self, prepared_request, response):
        host = requests.utils.urlparse(prepared_request.url).hostname
        if host != URS and "Authorization" in prepared_request.headers:
            del prepared_request.headers["Authorization"]


def fetch_all():
    s = EarthdataSession()
    out = []
    for ds, names in FILES.items():
        for n in names:
            out.append(download(f"{BASE}/{ds}/data/{n}", raw_dir("shift") / ds / n, session=s))
    return out
