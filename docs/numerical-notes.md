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

**What.** Pure-Python pyrobustlm is roughly 5-50x slower than R's C-backed
lmrob on default settings. Bench (single thread, AMD x86_64, OpenBLAS):

| n | p | R | pyrobustlm | ratio |
|---|---|---|------------|-------|
| 100  | 5  | 8 ms   | 390 ms | 49x |
| 500  | 10 | 30 ms  | 287 ms | 9.6x |
| 1000 | 10 | 46 ms  | 373 ms | 8.1x |
| 2000 | 20 | 249 ms | 574 ms | 2.3x |

**Why.** Python-level overhead in the resampling loop dominates on small
problems. On larger problems NumPy/BLAS catches up.

**Future work.** Phase 11 plans Cython acceleration of the inner loops
(``_core/_fast_s.pyx``). Plan target: ≤1.3x R on n>=1e4, multi-threaded
beat R on n>=1e3.

### 3. ``init="M-S"`` not yet implemented (Phase 5)

**What.** ``Control(init="M-S")`` raises ``NotImplementedError``. Designs
with categorical predictors that produce frequently singular subsamples
should bump ``Control(mts=...)`` higher (default 1000 already covers our
reference corpus). Phase 5 will implement Maronna-Yohai 2000 properly.

### 4. ``vcov_w`` falls back to ``vcov_avar1`` (Phase 7)

**What.** ``Control(cov=".vcov.w")`` currently emits a ``RuntimeWarning``
and uses the ``avar1`` sandwich. The KS2011 finite-sample corrections will
land in a follow-up.
