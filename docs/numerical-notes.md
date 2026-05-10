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

**What.** ``init="M-S"`` and ``init="auto"`` use a full port of
robustbase's ``R_lmrob_M_S`` in four phases (orthogonalize via L1,
subsample many candidates in orth space, transform back, descent via
alternating L1 + weighted-LS).

**Gap vs R.** On the ``education`` reference (Y ~ Region + X1 + X2 + X3)
at ``nResample=2000``::

    R   coef: (-135.7, -20.6, -9.9, 24.6, 0.034, 0.043, 0.579)  scale=26.4
    py  coef: (-138.2, -19.4, -10.0, 24.2, 0.035, 0.043, 0.584)  scale=26.4

Coefficient max-rerr ~1.8e-2, scale rerr ~6e-3. Most of the residual gap
is the usual RNG-basin drift (PCG64 vs MT) rather than an algorithmic
mismatch. The gap shrinks with more resamples; R's default of 500 will
land slightly different from ours.

**User guidance.** ``init="S"`` matches R to ``rtol=1e-6`` on continuous
designs. For factor designs, ``init="M-S"`` matches R within ~2% on
default tunings (use ``Control(nResample=2000)`` for tighter agreement).

### 4. ``vcov_w`` (Phase 7)

**What.** ``Control(cov=".vcov.w")`` now implements all five
``cov.corrfact`` branches (``asympt``, ``empirical``, ``tau``,
``hybrid``, ``tauold``), the five ``cov.dfcorr`` modes, the three
``cov.resid`` modes, and the Huber finite-sample correction. R's
setting-driven defaults are honoured (``"D" in method`` triggers
``cov.hubercorr=False`` and ``cov_corrfact="tau"``).

On stackloss with ``setting="KS2014"`` (and ``"KS2011"``), the cov
matrix matches R element-wise to ``rtol=1e-3``.

### 5. ``anova()`` chained mode

**What.** Calling ``anova(m1, m2, m3, ...)`` compares each adjacent pair
sequentially: ``m2`` vs ``m1``, then ``m3`` vs ``m2``, etc. R's
``anova.lmrob`` chained mode keeps the largest model as the reference for
every row, but its current implementation has a bug: row 3 onwards prints
the same statistic as row 2 because the iteration's reference state
collapses to the previous reduced fit's term set after the first pair.

**Why acceptable.** Pair-wise calls ``anova(full, reduced)`` match R
element-wise (chi-sq and p-value to ``rtol=2e-3`` on stackloss). Our
sequential chained behaviour matches base R's ``anova.lm`` semantics,
which is the conventional reading of nested-model anova tables.

**Where.** ``tests/validation/test_summary_anova.py::test_anova_*``.

### 6. ``vcov_avar1`` matches R element-wise

**What.** With the corrected ``Mchi(deriv=1) = chi'`` mapping (R's chi is
normalised so ``chi(inf) = 1``; ``chi' = (1/rho_unnorm(inf)) * psi``) and
the proper R formula ported from ``lmrob.MM.R:510-577``, the covariance
matrix now matches R element-wise to ``rtol=1e-3`` on stackloss/delivery/
phosphor.
