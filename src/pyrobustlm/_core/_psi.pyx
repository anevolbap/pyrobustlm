# cython: language_level=3
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Cython-accelerated bisquare psi/wgt/rho/psi_prime + the m_scale iteration.
# Used by the resampling loop in fast_s for the (extremely common) bisquare
# default. Pure-NumPy implementations in pyrobustlm._psifuns remain canonical.

cimport cython
from libc.math cimport fabs, sqrt
import numpy as np
cimport numpy as cnp

cnp.import_array()


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bisquare_psi(double[::1] x, double k, double[::1] out):
    """Tukey bisquare psi: x * (1 - (x/k)^2)^2 for |x| <= k, else 0."""
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double a, u
    for i in range(n):
        a = x[i] / k
        if a < -1.0 or a > 1.0:
            out[i] = 0.0
        else:
            u = 1.0 - a * a
            out[i] = x[i] * u * u


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bisquare_wgt(double[::1] x, double k, double[::1] out):
    """Bisquare weight: psi(x)/x = (1 - (x/k)^2)^2 for |x| <= k, else 0."""
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double a, u
    for i in range(n):
        a = x[i] / k
        if a < -1.0 or a > 1.0:
            out[i] = 0.0
        else:
            u = 1.0 - a * a
            out[i] = u * u


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bisquare_rho(double[::1] x, double k, double[::1] out):
    """Normalised chi: t * (3 - 3t + t^2) for t = (x/k)^2 <= 1, else 1."""
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double t, ax
    for i in range(n):
        ax = x[i] if x[i] >= 0 else -x[i]
        if ax > k:
            out[i] = 1.0
        else:
            t = (x[i] / k) * (x[i] / k)
            out[i] = t * (3.0 + t * (-3.0 + t))


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bisquare_psi_prime(double[::1] x, double k, double[::1] out):
    """psi'(x; k) = (1 - (x/k)^2)(1 - 5(x/k)^2)  for |x| <= k, else 0."""
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double a, a2
    for i in range(n):
        a = x[i] / k
        if a < -1.0 or a > 1.0:
            out[i] = 0.0
        else:
            a2 = a * a
            out[i] = (1.0 - a2) * (1.0 - 5.0 * a2)


# ---------------------------------------------------------------------------
# M-scale inner iteration for bisquare. Returns the converged scale.
# ---------------------------------------------------------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def m_scale_bisquare(
    double[::1] r,
    double k,
    double b0,
    double init_scale,
    int max_iter,
    double tol,
    int p,
):
    """One M-scale iteration loop for bisquare chi.

    Returns the converged sigma satisfying ``mean_{n-p} chi(r/sigma) = b0``.
    """
    cdef Py_ssize_t i, it, n = r.shape[0]
    cdef double s = init_scale
    cdef double prev = init_scale
    cdef double a, ax, t, sum_chi, diff
    cdef double inv_npmp = 1.0 / (n - p)
    if s <= 0.0:
        return 0.0
    for it in range(max_iter):
        sum_chi = 0.0
        for i in range(n):
            a = r[i] / s
            ax = a if a >= 0.0 else -a
            if ax > k:
                sum_chi += 1.0
            else:
                t = (a / k) * (a / k)
                sum_chi += t * (3.0 + t * (-3.0 + t))
        s = s * sqrt(sum_chi * inv_npmp / b0)
        diff = s - prev if s >= prev else prev - s
        if diff <= tol * prev:
            return s
        prev = s
    return s
