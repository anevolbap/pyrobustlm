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

**What.** With ``seed`` held fixed, our default fast-S
(``Control(rng="PCG64")``) converges to slightly different S
coefficients than R because we use NumPy's PCG64 instead of R's
Mersenne Twister. The set of resampled p-subsets differs, so the
best-of-best_r candidate is sometimes a different starting point.
After MM, the final beta agrees with R's MM beta to about
``rtol=5e-5`` on stackloss; on small-n datasets (pension, starsCYG)
the divergence can be ``rtol=1e-1``.

For redescending psi (hampel, ggw) on n=21 stackloss the basin drift
also leaks into the parametric cov: a few rweights differ from R's by
about 5e-3 (in 4 of 21 observations near the bisquare/ggw boundary),
and that propagates through ``X^T diag(rweights) X`` to a cov diagonal
rerr of about 0.4-0.5 on the most sensitive elements. Coef and scale
stay within ~1e-3 rerr of R; only the parametric vcov is sensitive on
this regime. See ``psi_hampel`` and ``psi_ggw`` rows in
``docs/bench-report.md``.

**Opt-in fix.** ``Control(rng="R")`` (v0.5.16+) drives the resample
loop through ``pylmrob.r_set_seed`` + ``r_sample_noreplace`` /
``r_subsample_nonsingular``, byte-identical to robustbase's
``unif_rand`` stream. End-to-end fits now agree with R's ``lmrob``
across the 10-dataset classical corpus to **rtol=1e-5** on
coefficients and ``~3.2e-6`` on scale. Forces ``n_workers=1`` and
``engine_c=False``; see [`rng-r-perf`](rng-r-perf.md) for wall-clock
costs.

The residual ``~3.2e-6`` scale floor is from accumulated LAPACK
``dgels`` ULP differences across the IRWLS solves; it's identical
to pylmrob's historic PCG64-path scale floor and is **not** caused
by RNG or MAD. Closing it would require matching R's exact BLAS
implementation; see ``plan.md`` §11.3 for details.

#### Investigation log: confirmed irreducible

Closing the gap was explicitly attempted in v0.5.19 dev cycle. None of
the following moved pylmrob closer to R on stackloss with
``Control(rng="R")``:

- ``OPENBLAS_NUM_THREADS=1`` + ``OMP_NUM_THREADS=1`` (set before
  numpy/scipy import). No change.
- Tightening ``rel_tol`` from ``1e-7`` to ``1e-14`` and ``max_it`` to
  200. The fit converges to a slightly *different* fixed point
  (intercept gap grows from 1.69e-5 to 2.04e-5). Two stationary
  points exist; pylmrob's BLAS lands on one, R's on the other.
- Forcing ``initial_scale=-1`` in ``cy_refine_fast_s_r`` so the
  survivor refinement re-MADs the residuals (matches R's call site
  exactly). Made the gap **worse** (3.0e-5 vs 1.7e-5).

The MM step is already routed through ``cy_lmrob_mm`` (LAPACK
``dgels``, the same QR-based routine ``robustbase::rwls`` uses) when
``rng="R"``; that ruled out gelsd/dgels divergence as the source. The
convergence test is also bit-identical to R's
``d_beta <= epsilon * fmax(epsilon, ||beta_new||_1)``.

**Conclusion.** The ``~1.7e-5`` intercept gap on stackloss is the
floor for pylmrob's BLAS environment. Users who need sub-1e-5 R
agreement should call R via rpy2; pylmrob's ``rng="R"`` is
production-ready for rtol=1e-5 reproducibility.

**Why acceptable for default.** Plan §5.2 documents the default RNG
strategy. The PCG64 cov drift is a function of how close the n=21
stackloss observations sit to the psi-redescending region, not an
algorithmic mismatch; the median cov diag rerr across the 34-case
corpus is 7.95e-07.

**Where.** Tests at ``tests/validation/test_vs_r_classical.py`` apply
per-dataset tolerances reflecting this sensitivity (PCG64 path).
``tests/validation/test_lmrob_rng_r_vs_R.py`` exercises the tighter
``rng="R"`` parity.

### 2. Performance vs R

**Default ``Control()``** (since v0.5.11) routes the entire fit
through a monolithic Cython kernel. Median wall-clock across the
34-case bench corpus is **0.93x R**; pylmrob is faster than R on
more than half of cases. Pass ``Control(engine_c=False)`` to force
the legacy NumPy path (median 2.74x R).

| dataset / setting | engine_c=False | default | ratio vs R |
|---|---|---|---|
| stackloss MM     | 25 ms | 4.7 ms | 1.2x R |
| stackloss KS2014 | 71 ms | 6.0 ms | 0.7x R |
| phosphor MM      | 26 ms | 4.3 ms | 1.4x R |
| salinity MM      | 28 ms | 4.2 ms | 0.8x R |
| delivery MM      | 27 ms | 4.2 ms | 1.0x R |

**Why this is the floor.** The Cython kernel performs the same BLAS
work R's C kernel does (we share OpenBLAS). The remaining gap is one
Python/C boundary cross at the top of every fit plus formulaic
parsing for non-trivial formulas. Going below R wall-clock would
mean dropping the Python public API.

**The engine_c trade-off** is byte-level RNG drift; see
``docs/engine_c.md`` for the long explanation.

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

### 5. Cython resampling kernel (bisquare hot path)

**What.** ``pylmrob._core._fast_s.cy_resample_iter_bisquare`` runs one
fast-S resampling iteration in ``nogil`` C: dgesv for the p-subset solve,
dgels for the IRWLS step, and an inlined m_scale iteration matching
``lmrob.c::find_scale``. The bisquare path of ``fast_s`` dispatches to
this kernel automatically.

**Speedup.** Serial vs the v0.2.0 NumPy implementation, all on a 16-core
OpenBLAS Linux box with BLAS pinned to 1 thread (so threading benefits
are clean):

| n / p / nResample | NumPy serial | Cython serial | Cython 8-thread |
|---|---|---|---|
| 100 / 5 / 500 | 102 ms | 43 ms | 42 ms |
| 500 / 10 / 500 | 261 ms | 128 ms | 68 ms |
| 1000 / 10 / 500 | 309 ms | 225 ms | 82 ms |
| 2000 / 20 / 500 | 634 ms | 562 ms | 187 ms |
| 5000 / 30 / 2000 | 9.7 s | 7.2 s | 2.0 s |

R wall-clock at n=2000/p=20 was 249 ms in the v0.1.0 benchmark; we
are now faster at that size with 8 threads.

**Where the rest of the gap lives.** At n=100 the iteration body is
small enough that subset-draw RNG and the survivor-refinement loop
(both still pure Python) dominate. Closing them would mean either
moving RNG into Cython (``np.random.cython`` extension or a struct-based
PCG64) or porting ``_refine_to_convergence`` to the same kernel. Both
deferred.

**Where.** ``src/pylmrob/_core/_fast_s.pyx`` and
``tests/unit/test_fast_s_parallel.py``.

### 6. Thread-based parallel resampling

**What.** ``Control(n_workers=...)`` opts the fast-S resampling loop into a
``ThreadPoolExecutor``. Each worker draws from a ``SeedSequence``-spawned
PCG64, so results are deterministic for a given
``(seed, n_workers, nResample)`` triple. ``n_workers=1`` (the default) is
serial and bit-identical with pre-parallel releases. ``n_workers=0`` is
auto, only enabling threading when ``n * p^2 >= 1e6`` and ``nResample >=
250`` (heuristic measured on a 16-core OpenBLAS Linux box).

**Why not OpenMP/Cython prange.** The resampling iteration is dominated
by NumPy's ``np.linalg.solve``, ``lstsq``, ``svd``, and matmul. Those
already release the GIL, so a thread pool is the realistic equivalent
without rewriting the loop body in C+LAPACK. True OpenMP would require
porting the linalg path off NumPy.

**Speedup numbers** (16-core, OpenBLAS pinned to 1 thread, n=5000, p=30,
nResample=2000): 9.74 s serial, 5.77 s auto. About 1.7x.

**Why small-n doesn't speed up.** For tiny problems the per-iteration
work is mostly Python (RNG draws, dictionary lookups, heap updates),
which holds the GIL. Threading there pays the pool overhead without
freeing real CPU. The 12x R gap at n=100 is GIL-bound, not BLAS-bound;
closing it would require pushing the loop body into Cython.

**Where.** ``tests/unit/test_fast_s_parallel.py``.

### 7. ``anova()`` chained mode

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

### 8. ``vcov_avar1`` matches R element-wise

**What.** With the corrected ``Mchi(deriv=1) = chi'`` mapping (R's chi is
normalised so ``chi(inf) = 1``; ``chi' = (1/rho_unnorm(inf)) * psi``) and
the proper R formula ported from ``lmrob.MM.R:510-577``, the covariance
matrix now matches R element-wise to ``rtol=1e-3`` on stackloss/delivery/
phosphor.

### 9. D-step ``kappa``: match R's quadrature, not the integral

**What.** ``robustbase:::lmrob.kappa`` solves
``E[psi(Z) Z - kappa wgt(Z)] = 0`` with ``uniroot``. The expression is
linear in ``kappa``, so the root is ``E[psi(Z) Z] / E[wgt(Z)]`` and the
only question is how the expectations are taken. R takes them with
``robustbase:::lmrob.E``, which does **not** integrate: it applies an
``numpoints``-node Gauss-Hermite rule, and ``lmrob.control`` sets
``numpoints = 10``.

Ten nodes resolve the kinked families only roughly, so R's ``kappa``
sits up to 8e-3 away from the exact integral:

| family | R (GH, 10 nodes) | exact integral | rerr |
|---|---|---|---|
| bisquare | 0.8280907302 | 0.8280771566 | 1.6e-05 |
| hampel | 0.8504151400 | 0.8569775806 | 7.7e-03 |
| optimal | 0.9361829182 | 0.9355077953 | 7.2e-04 |
| lqq | 0.8618074039 | 0.8626400360 | 9.7e-04 |
| ggw (b=1) | 0.8989804359 | 0.8914986546 | 8.3e-03 |
| ggw (b=1.5) | 0.8590698461 | 0.8597035165 | 7.4e-04 |

We previously computed the exact integral with ``scipy.integrate.quad``
and inherited the whole right-hand column as error.

**Why the rule and not the integral.** The goal is agreement with R.
R's D-scale is defined by the value R actually uses, so reproducing the
quadrature is reproducing the estimator. ``Control.numpoints`` exists so
the node count stays a parameter rather than a hidden constant, and
``tests/unit/test_d_scale_kappa.py`` pins both the R values and the fact
that the exact integral genuinely differs, so a later reader does not
"fix" it back.

**A second, larger error on the same line.** The Cython D-step carried
its own hardcoded ``kappa`` table, copied from the old ``quad`` output.
Its ggw case-4 entry was a copy of the case-1 value, so every ggw fit
through the default engine ran with a ``kappa`` 3.8% wrong. ``kappa`` is
now computed once per fit in Python and passed into the kernel; there is
no table left to drift, and non-default tunings work instead of silently
falling back.

End-to-end KS2014 scale on stackloss vs R, worst over five seeds:

| psi | before | after |
|---|---|---|
| bisquare | 5.0e-01 | 1.2e-07 |
| hampel | 7.7e-03 | 1.1e-07 |
| optimal | 7.2e-04 | 4.4e-12 |
| lqq | 9.7e-04 | 8.6e-08 |
| ggw | 3.8e-02 | 7.1e-07 |

The bisquare row was a different fault found in the same investigation;
see entry 10.

**Where.** ``pylmrob/d_scale.py::kappa``,
``tests/unit/test_d_scale_kappa.py``,
``tests/integration/test_kernel_parity.py::test_cy_d_scale_matches_numpy``.

### 10. Degenerate zero-scale candidates in the Cython fast-S

**What.** ``_mscale_generic`` can return a non-positive scale for a
resampled candidate. Zero is smaller than every real scale, so such a
candidate always won the best-of-``best_r`` comparison. ``cy_lmrob_fit``
then returned ``status = 0`` (success) with ``scale = 0`` and never wrote
``beta_init_out``, and the caller read that uninitialised buffer:
coefficients around 1e241 fed straight into ``vcov_avar1``.

It needed roughly 1000 resamples to hit, so it looked seed-dependent: on
stackloss with ``setting="KS2014"`` and ``psi="bisquare"`` it fired on 3
of 4 seeds, giving scales 36-50% away from R while the NumPy path was
correct on all four.

**Why these are not exact fits.** A zero M-scale is a legitimate S
solution when the fit passes exactly through more than half the data.
That is not the case here: over all 5719 non-singular 4-subsets of
stackloss the best any candidate achieves is 8 exactly-zero residuals,
and ``(n - p) * b0 = 8.5`` are required. Checked in double precision
from both plausible starting scales (MAD and ``max|r| / k0``), no subset
collapses. These are numerical degeneracies, so they are now rejected
and the search continues. A genuine exact fit is still detected earlier,
by ``max|r| == 0``.

**Knock-on.** The ``FloatingPointError`` retry in ``lmrob()`` existed to
paper over the singular ``X'WX`` this produced. Over a 500-fit sweep (10
classical datasets x 5 psi families x 10 seeds) it now fires zero times.
It is kept for pathological data but warns, because falling back changes
the estimator and that should not be silent.

**Where.** ``pylmrob/_core/_psi_kernels.pxi``,
``pylmrob/_core/_lmrob.pyx::cy_lmrob_fit``,
``tests/integration/test_engine_c_parity.py``.
