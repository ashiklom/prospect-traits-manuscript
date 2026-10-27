"""Leaf radiative transfer models.

Every model is implemented twice: in NumPy (the reference, used for tests and
for simulating spectra from posterior draws) and in PyTensor (used inside the
PyMC inversion). Origins and citations are listed in ``data/PROVENANCE.md``.
"""
