# cython: language_level=3
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Cython port of pylmrob.rng's hot path: R-compatible MT19937 tempering,
# r_sample_noreplace, r_subsample_nonsingular. Each function takes the
# state buffer (uint32[624]) and current position and updates them in
# place. All inner loops are nogil.

cimport cython
from libc.math cimport fabs
from libc.stdint cimport int64_t, uint32_t
from libc.stdlib cimport free, malloc

import numpy as np
cimport numpy as cnp

cnp.import_array()


cdef int _N = 624
cdef int _M = 397
cdef uint32_t _MATRIX_A = <uint32_t>0x9908B0DFu
cdef uint32_t _UPPER_MASK = <uint32_t>0x80000000u
cdef uint32_t _LOWER_MASK = <uint32_t>0x7FFFFFFFu
cdef uint32_t _TEMPER_B = <uint32_t>0x9D2C5680u
cdef uint32_t _TEMPER_C = <uint32_t>0xEFC60000u


cdef inline void _regenerate(uint32_t* s) nogil:
    cdef int i
    cdef uint32_t y
    for i in range(_N - _M):
        y = (s[i] & _UPPER_MASK) | (s[i + 1] & _LOWER_MASK)
        s[i] = s[i + _M] ^ (y >> 1) ^ (_MATRIX_A if (y & 1) else 0)
    for i in range(_N - _M, _N - 1):
        y = (s[i] & _UPPER_MASK) | (s[i + 1] & _LOWER_MASK)
        s[i] = s[i + _M - _N] ^ (y >> 1) ^ (_MATRIX_A if (y & 1) else 0)
    y = (s[_N - 1] & _UPPER_MASK) | (s[0] & _LOWER_MASK)
    s[_N - 1] = s[_M - 1] ^ (y >> 1) ^ (_MATRIX_A if (y & 1) else 0)


cdef inline uint32_t _next_uint32(uint32_t* s, int* pos) nogil:
    cdef uint32_t y
    if pos[0] >= _N:
        _regenerate(s)
        pos[0] = 0
    y = s[pos[0]]
    pos[0] += 1
    y ^= y >> 11
    y ^= (y << 7) & _TEMPER_B
    y ^= (y << 15) & _TEMPER_C
    y ^= y >> 18
    return y


cdef inline double _unif_rand(uint32_t* s, int* pos) nogil:
    cdef uint32_t y = _next_uint32(s, pos)
    cdef double u = <double>y * (1.0 / 4294967296.0)  # 2^-32
    if u <= 0.0:
        return 0.5 / 4294967296.0
    if u >= 1.0:
        return 1.0 - 0.5 / 4294967296.0
    return u


@cython.boundscheck(False)
@cython.wraparound(False)
def cy_r_unif_rand_n(uint32_t[::1] state, int pos, int n):
    """Pull ``n`` unif_rand draws nogil. Returns ``(out, new_pos)``."""
    cdef cnp.ndarray[double, ndim=1, mode="c"] out = np.empty(n, dtype=np.float64)
    cdef double* out_ptr = <double*>cnp.PyArray_DATA(out)
    cdef uint32_t* s_ptr = &state[0]
    cdef int p = pos
    cdef int i
    with nogil:
        for i in range(n):
            out_ptr[i] = _unif_rand(s_ptr, &p)
    return out, p


@cython.boundscheck(False)
@cython.wraparound(False)
def cy_r_sample_noreplace(uint32_t[::1] state, int pos, int n, int k):
    """Knuth swap-and-replace, byte-identical to robustbase's sample_noreplace.

    Returns ``(out_indices, new_pos)``. Mutates ``state`` in place when the
    MT regenerate step fires.
    """
    if k < 0 or k > n:
        raise ValueError("need 0 <= k <= n")
    cdef cnp.ndarray[int64_t, ndim=1, mode="c"] out = np.empty(k, dtype=np.int64)
    cdef int64_t* out_ptr = <int64_t*>cnp.PyArray_DATA(out)
    cdef int64_t* ind_space = <int64_t*>malloc(n * sizeof(int64_t))
    if ind_space == NULL:
        raise MemoryError("cy_r_sample_noreplace: out of memory")
    cdef uint32_t* s_ptr = &state[0]
    cdef int p = pos
    cdef int i, j, nn
    cdef double u
    with nogil:
        for i in range(n):
            ind_space[i] = i
        nn = n
        for i in range(k):
            u = _unif_rand(s_ptr, &p)
            j = <int>(nn * u)
            out_ptr[i] = ind_space[j]
            nn -= 1
            ind_space[j] = ind_space[nn]
    free(ind_space)
    return out, p


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cy_r_subsample_nonsingular(
    uint32_t[::1] state,
    int pos,
    double[:, ::1] X,
    int p,
    int mts,
    double tol_inv,
):
    """Pick a non-singular p-row subset of X, byte-identical to robustbase's
    subsample() with ss=1. Returns ``(idc_or_None, new_pos)``.
    """
    cdef int n = X.shape[0]
    cdef int m = X.shape[1]
    if p < 1 or p > m:
        raise ValueError("need 1 <= p <= ncol(X)")
    if p > n:
        raise ValueError("need p <= nrow(X)")

    cdef cnp.ndarray[int64_t, ndim=1, mode="c"] idc_arr = np.empty(p, dtype=np.int64)
    cdef int64_t* idc = <int64_t*>cnp.PyArray_DATA(idc_arr)
    cdef int64_t* perm = <int64_t*>malloc(n * sizeof(int64_t))
    cdef int64_t* ind_space = <int64_t*>malloc(n * sizeof(int64_t))
    cdef int64_t* idr = <int64_t*>malloc(p * sizeof(int64_t))
    cdef double* lu = <double*>malloc(p * p * sizeof(double))
    cdef double* v = <double*>malloc(p * sizeof(double))
    if perm == NULL or ind_space == NULL or idr == NULL or lu == NULL or v == NULL:
        free(perm); free(ind_space); free(idr); free(lu); free(v)
        raise MemoryError("cy_r_subsample_nonsingular: out of memory")

    cdef uint32_t* s_ptr = &state[0]
    cdef double* X_ptr = &X[0, 0]
    cdef int pp = pos
    cdef int attempt = 0
    cdef int i, j, k, kk, ll, mu, nn
    cdef int sing, failed, success, found
    cdef double tmpd, ss, pivot, tmp_swap, u
    cdef int64_t tmpi

    found = 0

    with nogil:
        while True:
            # Step 1: full permutation via sample_noreplace into perm[].
            for i in range(n):
                ind_space[i] = i
            nn = n
            for i in range(n):
                u = _unif_rand(s_ptr, &pp)
                j = <int>(nn * u)
                perm[i] = ind_space[j]
                nn -= 1
                ind_space[j] = ind_space[nn]

            # Step 2: incremental LU with column pivot.
            for k in range(p):
                idr[k] = k
            i = 0
            failed = 0
            success = 1
            for j in range(p):
                sing = 1
                while sing:
                    if i + j >= n:
                        failed = 1
                        break
                    idc[j] = perm[i + j]
                    if j == 0:
                        for k in range(p):
                            v[k] = X_ptr[idc[j] * m + idr[k]]
                    else:
                        for k in range(j):
                            lu[k * p + j] = X_ptr[idc[j] * m + idr[k]]
                        for kk in range(j):
                            ss = lu[kk * p + j]
                            for ll in range(kk):
                                ss -= lu[kk * p + ll] * lu[ll * p + j]
                            lu[kk * p + j] = ss
                        for k in range(j, p):
                            ss = X_ptr[idc[j] * m + idr[k]]
                            for ll in range(j):
                                ss -= lu[k * p + ll] * lu[ll * p + j]
                            v[k] = ss
                    if j < p - 1:
                        tmpd = fabs(v[j])
                        mu = j
                        for k in range(j + 1, p):
                            if fabs(v[k]) > tmpd:
                                mu = k
                                tmpd = fabs(v[k])
                        if tmpd >= tol_inv:
                            tmp_swap = v[j]; v[j] = v[mu]; v[mu] = tmp_swap
                            tmpi = idr[j]; idr[j] = idr[mu]; idr[mu] = tmpi
                            pivot = v[j]
                            for k in range(j + 1, p):
                                lu[k * p + j] = v[k] / pivot
                            if j > 0 and mu != j:
                                for k in range(j):
                                    tmp_swap = lu[j * p + k]
                                    lu[j * p + k] = lu[mu * p + k]
                                    lu[mu * p + k] = tmp_swap
                    if fabs(v[j]) < tol_inv:
                        i += 1
                        continue
                    lu[j * p + j] = v[j]
                    sing = 0
                if failed:
                    success = 0
                    break

            if success:
                found = 1
                break
            attempt += 1
            if attempt >= mts:
                found = 0
                break

    free(perm); free(ind_space); free(idr); free(lu); free(v)
    if found:
        return idc_arr, pp
    return None, pp
