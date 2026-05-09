# Numerical notes

Documented divergences from R's `robustbase::lmrob`.

When a Python output deviates from R beyond the agreed tolerance
(`tests/conftest.py::DEFAULT_TOLERANCES`), the choice is to either
(a) fix the bug or (b) document the divergence here with rationale.

For each entry record:

- **What** differs (which output, by how much, on what input).
- **Why** the divergence is acceptable (e.g. unavoidable RNG drift, a
  numerically more stable algorithm, an upstream R bug we chose not to
  reproduce).
- **Where** to find the comparing test and the relevant source line.

---

## Entries

### 1. Initial-S basin sensitivity to RNG (Phase 4)

**What.** With ``seed`` held fixed, our fast-S converges to slightly different
S coefficients than R because we use NumPy's PCG64 instead of R's Mersenne
Twister. The set of resampled p-subsets differs, so the best-of-best_r
candidate is sometimes a different starting point. After MM, the final
beta agrees with R's MM beta to about ``rtol=5e-5`` on stackloss; on
small-n datasets (pension, starsCYG) the divergence can be ``rtol=1e-1``.

**Why acceptable.** Plan §5.2 explicitly waives bit-identical RNG
reproducibility with R.

**Where.** Tests at ``tests/validation/test_vs_r_classical.py`` apply
per-dataset tolerances reflecting this sensitivity.

### 2. Performance vs R (Phase 11)

**What.** With the bisquare Cython kernel for psi/wgt/rho/m_scale wired into
the hot path of ``fast_s`` and ``m_scale``, pyrobustlm is now within 1.1x to
12x of R's C-backed lmrob on the bisquare-default path:

| n | p | R | pyrobustlm | ratio |
|---|---|---|------------|-------|
| 100  | 5  | 8 ms   | 94 ms  | 12x  |
| 500  | 10 | 30 ms  | 79 ms  | 2.6x |
| 1000 | 10 | 46 ms  | 119 ms | 2.6x |
| 2000 | 20 | 249 ms | 277 ms | 1.1x |

**Why.** Cython kernels remove the Python-loop overhead in the inner
resampling loop. On large problems NumPy/BLAS dominates and we approach R.
On small problems Python's ``np.linalg.lstsq`` and ``np.linalg.solve``
overhead dominates.

**Remaining gap.** Other psi families (huber, hampel, optimal, lqq, ggw)
still use the NumPy path; Cython kernels for them follow the same pattern.

### 3. M-S estimator (Phase 5)

**What.** ``init="M-S"`` and ``init="auto"`` both work. The estimator is
a simplified Maronna-Yohai 2000: alternating L1 (via SciPy linprog) on the
factor block and S (via fast-S) on the continuous block.

**Gap vs R.** On the ``education`` reference (Y ~ Region + X1 + X2 + X3),
our simplified M-S converges to a different solution than R's
``robustbase::lmrob(init="M-S")``::

    R   coef: (-135.7, -20.6, -9.9, 24.6, 0.034, 0.043, 0.579)  scale=26.4
    py  coef: (-158.6,  -8.2, -11.3, 20.1, 0.045, 0.043, 0.632)  scale=37.8

This is **not** RNG drift; it is an algorithmic gap. R's M-S
(``robustbase/src/lmrob.c::m_s_subsample`` + ``m_s_descent``, ~600 lines
of C) does many random subsample restarts followed by a careful descent
phase. We currently do a single L1/S alternating descent. Porting the
multi-restart + descent logic is tracked but not yet done.

**User guidance.** If your design is continuous, use ``init="S"`` (the
default) and you will match R within RNG drift. If you must use M-S,
treat the result as approximate and re-fit in R for publication-grade
inference.

### 4. ``vcov_w`` (Phase 7)

**What.** ``Control(cov=".vcov.w")`` implements the
``cov.corrfact = "asympt"`` branch of robustbase's ``.vcov.w`` (the
asymptotic-correction-factor version). The other branches (``empirical``,
``hybrid``, ``tau``, ``tauold``) and the optional Huber finite-sample
correction are deferred; they primarily affect inference in small
samples and do not change point estimates.

### 5. ``vcov_avar1`` matches R element-wise

**What.** With the corrected ``Mchi(deriv=1) = chi'`` mapping (R's chi is
normalised so ``chi(inf) = 1``; ``chi' = (1/rho_unnorm(inf)) * psi``) and
the proper R formula ported from ``lmrob.MM.R:510-577``, the covariance
matrix now matches R element-wise to ``rtol=1e-3`` on stackloss/delivery/
phosphor.
