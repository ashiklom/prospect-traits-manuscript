"""Common interface for leaf RTMs."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class RTM:
    """A leaf RTM evaluated on its native wavelength grid.

    Subclasses implement ``simulate`` (NumPy reference) and ``simulate_pt``
    (PyTensor); both take a mapping from parameter name to value and return
    ``(R, T)`` on ``self.wavelength``. ``params`` lists the parameters in the
    order used for inference; ``fixed`` holds values that are inputs to the
    RTM but are not inferred.
    """

    name: str
    description: str
    params: tuple[str, ...]
    wavelength: np.ndarray
    fixed: dict[str, float] = field(default_factory=dict)

    def theta(self, values) -> dict:
        """Merge ``values`` (dict, or sequence ordered like ``params``) with
        ``fixed``."""
        if not isinstance(values, dict):
            values = dict(zip(self.params, values, strict=True))
        missing = set(self.params) - set(values)
        if missing:
            raise ValueError(f"{self.name}: missing parameters {sorted(missing)}")
        return {**self.fixed, **values}

    def simulate(self, values, **kwargs):  # pragma: no cover - abstract
        raise NotImplementedError

    def simulate_pt(self, values, **kwargs):  # pragma: no cover - abstract
        raise NotImplementedError
