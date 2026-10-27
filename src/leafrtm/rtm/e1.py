"""Exponential integral E1 for the PROSPECT plate model.

Vendored from ashiklom/prospect-inversion-comparison @ 84504ed
(``python/prospect_np.py`` and ``python/prospect_pt.py``). Changes from the
original: the Swamee-Ojha approximation and the JAX kernels are dropped
(they were never the fastest option in that repo's benchmarks).

Two modes:

* ``exact``: Cephes ``exp1`` via SciPy.
* ``fast``: the NAG S13AAF Chebyshev expansions (Horner form) used in the
  original PROSPECT Fortran (PEcAnRTM ``prospect.expint.f90``). Max relative
  error vs. Cephes ~2e-12; changes PROSPECT R/T by < 6e-15.

The PyTensor Ops supply the analytic derivative dE1/dx = -exp(-x)/x and ship
native numba kernels so PyTensor's numba backend (and nutpie) never fall back
to object mode.
"""

from __future__ import annotations

import ctypes

import numba
import numpy as np
import pytensor.tensor as pt
from numba.extending import get_cython_function_address
from pytensor.graph.basic import Apply
from pytensor.graph.op import Op
from pytensor.link.numba.dispatch.basic import register_funcify_default_op_cache_key
from scipy.special import exp1

# NAG S13AAF Chebyshev expansions (Horner form).
E1_FAST_LO = (  # 0 < x <= 4, in t = x/2 - 1;  E1 = poly(t) - log(x)
    -3.60311230482612224e-13, 3.46348526554087424e-12, -2.99627399604128973e-11,
    2.57747807106988589e-10, -2.09330568435488303e-9, 1.59501329936987818e-8,
    -1.13717900285428895e-7, 7.55292885309152956e-7, -4.64980751480619431e-6,
    2.63830365675408129e-5, -1.37089870978830576e-4, 6.47686503728103400e-4,
    -2.76060141343627983e-3, 1.05306034687449505e-2, -3.57191348753631956e-2,
    1.07774527938978692e-1, -2.96997075145080963e-1, 8.64664716763387311e-1,
    7.42047691268006429e-1)
E1_FAST_HI = (  # x > 4, in t = 14.5/(x + 3.25) - 1;  E1 = exp(-x) * poly(t) / x
    -1.62806570868460749e-12, -8.95400579318284288e-13, -4.08352702838151578e-12,
    -1.45132988248537498e-11, -8.35086918940757852e-11, -2.13638678953766289e-10,
    -1.10302431467069770e-9, -3.67128915633455484e-9, -1.66980544304104726e-8,
    -6.11774386401295125e-8, -2.70306163610271497e-7, -1.05565006992891261e-6,
    -4.72090467203711484e-6, -1.95076375089955937e-5, -9.16450482931221453e-5,
    -4.05892130452128677e-4, -2.14213055000334718e-3, -1.06374875116569657e-2,
    -8.50699154984571871e-2, 9.23755307807784058e-1)


def _horner(coefs, t):
    y = np.zeros_like(t)
    for c in coefs:
        y = y * t + c
    return y


def e1_fast_np(x):
    """Fast, accurate E1(x), x > 0 (NAG S13AAF Chebyshev expansions)."""
    x = np.asarray(x, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        lo = _horner(E1_FAST_LO, 0.5 * x - 1.0) - np.log(x)
        hi = np.exp(-x) * _horner(E1_FAST_HI, 14.5 / (x + 3.25) - 1.0) / x
    return np.where(x <= 4.0, lo, hi)


E1_NP = {"exact": exp1, "fast": e1_fast_np}


class _E1Base(Op):
    """Elementwise E1(x) with the analytic gradient dE1/dx = -exp(-x)/x."""

    __props__ = ()

    def make_node(self, x):
        x = pt.as_tensor_variable(x)
        return Apply(self, [x], [x.type()])

    def infer_shape(self, node, shapes):
        return shapes

    def pullback(self, inputs, outputs, cotangents):
        (x,) = inputs
        (gz,) = cotangents
        return [gz * (-pt.exp(-x) / x)]


class E1Exact(_E1Base):
    """Exact E1 (Cephes)."""

    def perform(self, node, inputs, output_storage):
        output_storage[0][0] = np.asarray(exp1(inputs[0]), dtype=node.outputs[0].dtype)


class E1Fast(_E1Base):
    """Fast E1 (NAG S13AAF Chebyshev expansions)."""

    def perform(self, node, inputs, output_storage):
        output_storage[0][0] = np.asarray(e1_fast_np(inputs[0]), dtype=node.outputs[0].dtype)


e1_exact = E1Exact()
e1_fast = E1Fast()
E1_PT = {"exact": e1_exact, "fast": e1_fast}

_exp1_c = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)(
    get_cython_function_address("scipy.special.cython_special", "__pyx_fuse_1exp1")
)


@register_funcify_default_op_cache_key(E1Exact)
def _numba_funcify_E1Exact(op, node=None, **kwargs):
    @numba.njit(error_model="numpy")
    def e1_numba(x):
        out = np.empty_like(x)
        xf = x.ravel()
        of = out.ravel()
        for i in range(xf.size):
            of[i] = _exp1_c(xf[i])
        return out

    # The integer salts PyTensor's on-disk numba cache key; bump it whenever
    # the kernel changes (the default key ignores the kernel's code).
    return e1_numba, 2


_LO, _HI = E1_FAST_LO, E1_FAST_HI


@register_funcify_default_op_cache_key(E1Fast)
def _numba_funcify_E1Fast(op, node=None, **kwargs):
    # error_model="numpy": x = 0 or inf gives inf/nan (-> rejected NUTS
    # proposal) instead of raising ZeroDivisionError; no fastmath, since it
    # would let LLVM assume inputs are finite.
    @numba.njit(error_model="numpy")
    def e1_fast_numba(x):
        out = np.empty_like(x)
        xf = x.ravel()
        of = out.ravel()
        for i in range(xf.size):
            xi = xf[i]
            if xi <= 4.0:
                t = 0.5 * xi - 1.0
                y = 0.0
                for c in _LO:
                    y = y * t + c
                of[i] = y - np.log(xi)
            else:
                t = 14.5 / (xi + 3.25) - 1.0
                y = 0.0
                for c in _HI:
                    y = y * t + c
                of[i] = np.exp(-xi) * y / xi
        return out

    return e1_fast_numba, 2
