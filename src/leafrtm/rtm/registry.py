"""Registry of available leaf RTMs (constructed lazily and cached)."""

from __future__ import annotations

from functools import cache

from leafrtm.rtm.base import RTM
from leafrtm.rtm.leafsip import LeafSIP
from leafrtm.rtm.liberty import Liberty, LibertyFixed
from leafrtm.rtm.prospect import VERSIONS, Prospect

MODELS = (*VERSIONS, "liberty", "liberty_fixed", "leafsip")


@cache
def get_model(name: str) -> RTM:
    if name in VERSIONS:
        return Prospect.from_version(name)
    if name == "liberty":
        return Liberty.create()
    if name == "liberty_fixed":
        return LibertyFixed.create()
    if name == "leafsip":
        return LeafSIP.create()
    raise KeyError(f"Unknown RTM {name!r}; choose from {MODELS}")
