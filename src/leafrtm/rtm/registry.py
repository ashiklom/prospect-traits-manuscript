"""Registry of available leaf RTMs (constructed lazily and cached)."""

from __future__ import annotations

from functools import cache

from leafrtm.rtm.base import RTM
from leafrtm.rtm.liberty import Liberty
from leafrtm.rtm.prospect import VERSIONS, Prospect

MODELS = (*VERSIONS, "liberty")


@cache
def get_model(name: str) -> RTM:
    if name in VERSIONS:
        return Prospect.from_version(name)
    if name == "liberty":
        return Liberty.create()
    raise KeyError(f"Unknown RTM {name!r}; choose from {MODELS}")
