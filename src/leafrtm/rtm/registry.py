"""Registry of available leaf RTMs (constructed lazily and cached)."""

from __future__ import annotations

from functools import cache

from leafrtm.rtm.base import RTM
from leafrtm.rtm.prospect import VERSIONS, Prospect

MODELS = tuple(VERSIONS)


@cache
def get_model(name: str) -> RTM:
    if name in VERSIONS:
        return Prospect.from_version(name)
    raise KeyError(f"Unknown RTM {name!r}; choose from {MODELS}")
