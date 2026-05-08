# cython: language_level=3
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Phase 0 toolchain probe.
# Real kernels (_psi, _scale, _fast_s, _mm) replace this in their phases.

cimport cython
from libc.math cimport sqrt


def hello() -> str:
    """Return a string, proving the extension imports."""
    return "pyrobustlm._core._stub OK"


@cython.boundscheck(False)
@cython.wraparound(False)
def vec_norm(double[::1] x) -> float:
    """Euclidean norm of a contiguous double array.

    Smoke-test for typed memoryview compilation.
    """
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double s = 0.0
    for i in range(n):
        s += x[i] * x[i]
    return sqrt(s)
