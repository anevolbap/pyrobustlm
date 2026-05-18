# `pylmrob` — Python port of R `robustbase::lmrob`

> **Audience:** Claude Code (or any contributor) implementing this project end-to-end.
> **Goal:** A production-grade Python implementation of the `lmrob` MM-estimator with feature parity, numerically equivalent results (to converged tolerance), and runtime ≥ R's C-backed implementation.

---

## 0. North Star & Non-Goals

### 0.1 Goals (in priority order)

1. **Correctness first.** For converged outputs (final β, σ, weights, cov), match R's `robustbase::lmrob` to within `1e-6` relative tolerance on a fixed validation corpus.
2. **Feature parity.** Support all settings R supports: `setting ∈ {"KS2011", "KS2014", default MM}`, `init ∈ {"S", "M-S", "L1"}`, all six psi families (`bisquare`, `huber`, `hampel`, `optimal`, `ggw`, `lqq`), all covariance estimators (`.vcov.avar1`, `.vcov.w`, KS2014).
3. **Speed parity.** ≥ R on n ≥ 1e4; ≤ 1.5× R wall-clock on small problems. Beat R when multiple cores are available.
4. **Pythonic API.** `formulaic` formulas, NumPy/pandas inputs, scikit-learn-style estimator class, statsmodels-style results object.
5. **Clean packaging.** Cross-platform wheels (Linux/macOS/Windows × x86_64/arm64 × CPython 3.10–3.13).

### 0.2 Non-Goals (v1)

- Bit-identical reproducibility with R's RNG. We document this clearly.
- `lmrob..M..fit` standalone fitter, mixed-effects, GLM extensions, `nlrob`. (Stretch goals only.)
- Bayesian / posterior interfaces.

### 0.3 Licensing

R `robustbase` is **GPL-2 | GPL-3**. Two strategies:

- **Strategy A (recommended):** Reimplement from papers + R source as reference. License this project as **GPL-3** with prominent attribution to robustbase authors. Easiest legally, and aligns with downstream user expectations.
- **Strategy B:** Clean-room port from papers only, no reading of C source. Allows MIT/BSD licensing but is much slower and harder; strongly not recommended.

**Decision: Strategy A. License = GPL-3.** Add `LICENSE`, `NOTICE`, and per-file SPDX headers in Phase 0.

---

## 1. Required Reading (Do This First)

Claude Code: **before writing any code**, read and take notes on the following. Save notes to `docs/research-notes.md`.

### 1.1 Papers (in this order)

1. Yohai (1987) — "High breakdown-point and high efficiency robust estimates for regression". *Annals of Statistics* 15(2). Defines MM-estimators.
2. Salibian-Barrera & Yohai (2006) — "A fast algorithm for S-regression estimates". *JCGS* 15(2). The fast-S algorithm we will port.
3. Maronna & Yohai (2000) — "Robust regression with both continuous and categorical predictors". *J. Stat. Plan. Inf.* The M-S estimator for designs with factors.
4. Koller & Stahel (2011) — "Sharpening Wald-type inference in robust regression for small samples". *CSDA*. Defines the KS2011 setting and `.vcov.w`.
5. Koller & Stahel (2017) — "Nonsingular subsampling for regression S estimators with categorical predictors". *Comp. Stat.* Defines KS2014 default and `lqq` psi.
6. Koller (2016) — robustbase vignette `lmrob_simulation.pdf`. Reference benchmark numbers.

### 1.2 Source Code

Read in this order, with the goal of producing `docs/r-source-map.md` mapping every R function and C function to a planned Python module.

| Upstream (CRAN mirror `cran/robustbase`) | What it is | Read for |
|---|---|---|
| `R/lmrob.R` | Top-level `lmrob()` and dispatch | API surface |
| `R/lmrob.MM.R` | MM iteration | Phase 6 |
| `R/lmrob.M.S.R` | M-S initial | Phase 5 |
| `R/lmrobMMcontrol.R` | `lmrob.control` defaults | Phase 8 |
| `R/lmrobMMresid.R` | residuals, predict | Phase 8 |
| `src/lmrob.h` | constants, structs | Phase 2 onward |
| `src/lmrob.c` | `fast_s`, `fast_s_large_n`, `rwls`, `rho`, `psi` | Phase 4–6 |
| `src/lmrob-psifuns.c` | psi/chi/wgt/Epsi for all families | Phase 2 |
| `src/mc.c` | Mahalanobis-like utilities | reference only |
| `src/init.c` | `.Call` registrations | reference only |

### 1.3 Existing Python attempts (study, don't depend on)

- `https://github.com/deepak7376/robustbase` — only implements Sn/Qn/MAD/IQR scale estimators. Useful for those primitives' style; does **not** include `lmrob`.
- `robustbase-py-optimized` on PyPI — thin wrappers, incomplete, slow. Reference only.
- `statsmodels.robust.robust_linear_model.RLM` — M-estimator only, not MM. Not a substitute.

---

## 2. Tech Stack & Architecture Decisions

| Concern | Decision | Rationale |
|---|---|---|
| Hot kernels | **Cython** with typed memoryviews | True C speed, OpenMP `prange`, no JIT warm-up, clean wheel story |
| LAPACK access | `scipy.linalg.cython_lapack` (`dgels`, `dgeqrf`, `dgesv`, `dpotrf`) | Avoid Python overhead inside resampling loop |
| Parallelism | OpenMP via Cython `prange` on candidate loop | Matches R's `_OPENMP` parallelism |
| RNG | `numpy.random.Generator` (PCG64) per-thread, seeded from `SeedSequence` | Reproducible across runs; do not attempt R-MT bit-identity |
| Formulas | `formulaic` (fallback to `patsy`) | Modern, faster, actively maintained |
| Build | `meson-python` + Cython | Modern PEP 517 build for compiled wheels |
| Wheels | `cibuildwheel` in CI | Standard practice |
| Validation | `rpy2` driving R inside pytest | Run R and Python on identical data, compare |
| Numerical lib | NumPy + SciPy LAPACK | Sufficient and standard |
| Python target | 3.10+ | Matches NumPy / SciPy support windows |

### 2.1 Why not Numba / pure NumPy?

- **Numba:** Cold-start JIT compile time, weaker reproducibility, awkward to call BLAS from `@njit`, harder to ship as a binary wheel.
- **Pure NumPy:** Vectorization can't accelerate the inner resampling loop because each iteration solves a different small `p × p` system. R's C loop will dominate.

### 2.2 Package name

`pylmrob`. Reasoning: `robustbase` and `robustbase-py-optimized` are taken on PyPI; `pyrobust` is too generic. Pin name in `pyproject.toml` Phase 0.

---

## 3. Repository Layout

```
pylmrob/
├── pyproject.toml              # meson-python build, deps, metadata
├── meson.build                 # top-level build
├── LICENSE                     # GPL-3.0
├── NOTICE                      # robustbase attribution
├── README.md
├── PLAN.md                     # this file
├── CHANGELOG.md
├── .github/workflows/
│   ├── ci.yml                  # lint, type, test on Linux/macOS/Win
│   ├── wheels.yml              # cibuildwheel on tag
│   └── docs.yml                # Sphinx → GH Pages
├── docs/
│   ├── research-notes.md       # Phase 1 deliverable
│   ├── r-source-map.md         # mapping R/C functions → Python modules
│   ├── api/                    # Sphinx
│   └── numerical-notes.md      # known divergences from R, tolerances
├── src/pylmrob/
│   ├── __init__.py
│   ├── _version.py
│   ├── control.py              # Control dataclass, presets
│   ├── psi.py                  # Python wrappers for psi/chi/wgt
│   ├── scale.py                # M-scale wrappers
│   ├── lmrob.py                # public lmrob() + LmRob class
│   ├── results.py              # LmRobResults, summary table
│   ├── inference.py            # covariance estimators, Wald, CIs
│   ├── diagnostics.py          # plot(), influence, leverage
│   ├── formula.py              # formulaic glue
│   ├── ms_estimator.py         # M-S init for factor designs
│   ├── _utils.py               # design-matrix prep, rank checks
│   └── _core/                  # Cython
│       ├── __init__.py
│       ├── meson.build
│       ├── _psi.pyx            # psi/chi/wgt/Epsi/EDpsi tight loops
│       ├── _scale.pyx          # M-scale iteration
│       ├── _fast_s.pyx         # fast-S resampling + IRWLS
│       ├── _mm.pyx             # MM iteration
│       ├── _lapack_helpers.pyx # cython_lapack wrappers
│       └── _ctables.pyx        # ggw/lqq tuning constant tables
├── tests/
│   ├── conftest.py             # shared fixtures, R bridge
│   ├── data/                   # static datasets (stackloss, hbk, ...)
│   ├── reference/              # serialized R outputs (committed)
│   ├── unit/
│   │   ├── test_psi.py
│   │   ├── test_scale.py
│   │   ├── test_fast_s.py
│   │   ├── test_mm.py
│   │   ├── test_inference.py
│   │   └── test_control.py
│   ├── integration/
│   │   ├── test_lmrob_api.py
│   │   ├── test_factor_designs.py
│   │   └── test_predict.py
│   ├── validation/
│   │   ├── test_vs_r_classical.py
│   │   ├── test_vs_r_synthetic.py
│   │   └── test_vs_r_settings.py
│   ├── property/
│   │   └── test_equivariance.py # Hypothesis: affine equivariance
│   └── benchmarks/
│       ├── bench_fast_s.py
│       └── bench_full_fit.py
└── scripts/
    ├── generate_r_reference.R   # produces tests/reference/*.json
    ├── psi_constant_tables.R    # generate ggw/lqq tables → CSV
    └── profile.py
```

---

## 4. Phased Plan

Each phase has: **Inputs**, **Tasks** (each a single Claude Code session worth), **Acceptance**, **Validation Commands**, **Dependencies**.

---

### Phase 0 — Project Bootstrap (1–2 days)

**Inputs:** This document.

**Tasks:**

1. Initialize repo. Add `pyproject.toml` (meson-python), `meson.build`, `LICENSE` (GPL-3), `NOTICE` crediting robustbase authors (Maechler, Rousseeuw, Croux, Todorov, Ruckstuhl, Salibian-Barrera, Verbeke, Koller, Conceicao, Anna di Palma).
2. Configure tooling: `ruff` (lint+format), `mypy` (strict), `pytest`, `pytest-cov`, `pytest-xdist`, `hypothesis`, `pre-commit`.
3. Set up CI matrix (`.github/workflows/ci.yml`): Linux/macOS/Windows × Python 3.10/3.11/3.12/3.13. R + rpy2 install step. Cache R packages.
4. Add `cibuildwheel` workflow (triggered on tag).
5. Write `tests/conftest.py` with shared fixtures: `r_session` (rpy2), `rng` (deterministic), `tol` (numerical tolerances dict).
6. Write `scripts/generate_r_reference.R` skeleton.

**Acceptance:**
- `pip install -e .` succeeds locally.
- `pytest` runs (zero tests OK).
- CI is green.
- `import pylmrob; pylmrob.__version__` works.

**Validation:** `make ci-local` runs lint + type + test.

---

### Phase 1 — Reference Harness via R (2–3 days)

This is the testbed against which **everything** will be validated. Build it before any numerical code.

**Tasks:**

1. **`scripts/generate_r_reference.R`** — given a list of `(dataset_name, formula, control_args, seed)` tuples, runs `lmrob` in R and writes JSON containing: `coefficients`, `scale`, `weights`, `residuals`, `fitted.values`, `cov`, `df.residual`, `init.S$coefficients`, `init.S$scale`, `converged`, plus all psi/chi/wgt evaluations needed for unit tests.
2. **Datasets to include** (all bundled in `robustbase` itself, so easy to load): `stackloss`, `coleman`, `salinity`, `wood`, `hbk`, `starsCYG`, `delivery`, `education`, `aircraft`, `pension`, `phosphor`. Plus 6 synthetic configurations (n ∈ {50, 200, 1000}, p ∈ {3, 10}, contamination ∈ {0%, 10%, 30%}).
3. **Settings to include:** default MM, `setting="KS2011"`, `setting="KS2014"`, plus per-psi runs (`psi="bisquare"`, `"optimal"`, `"ggw"`, `"lqq"`).
4. **rpy2 bridge in `conftest.py`:** helper that, given the same `(dataset, formula, control)`, can either load the cached JSON reference or invoke R live.
5. **Tolerance config in `conftest.py`:** dict mapping output name → (rtol, atol). Start strict (1e-10 for psi, 1e-6 for converged β); loosen as needed with rationale logged.

**Acceptance:**
- `Rscript scripts/generate_r_reference.R` produces ≥ 50 reference JSON files, all loadable from Python.
- Loading helpers tested.

**Validation:**
```bash
Rscript scripts/generate_r_reference.R
pytest tests/conftest.py -k "test_reference_loader"
```

---

### Phase 2 — Psi / Chi / Weight Functions (1 week)

**Inputs:** Phase 1 reference values.

**The six families:** `huber`, `bisquare` (Tukey), `hampel`, `optimal`, `ggw` (generalized Gauss-weight), `lqq` (linear-quadratic-quadratic).

**For each family, six functions:** `rho(x, k)`, `psi(x, k)` = `rho'`, `psi_prime(x, k)` = `rho''`, `wgt(x, k)` = `psi(x)/x`, `Epsi2(k)` = E[ψ²], `EDpsi(k)` = E[ψ'].

**Tasks:**

1. **`_core/_psi.pyx`** — Cython implementations with `nogil` and `cpdef inline` scalar versions plus vectorized memoryview versions. One function per (family, kind). Use `libc.math` for `exp`, `fabs`, etc.
2. **For `ggw` and `lqq`:** these have no closed-form tuning constants. Port the precomputed tables from R (`R/lmrobMMcontrol.R` → `.psi.const`). Generate via `scripts/psi_constant_tables.R` and embed as `static const double[]` in `_ctables.pyx`.
3. **`psi.py`** — Python-facing wrappers: `psi(x, family, k)` dispatches to the Cython kernel. Handle scalar vs array inputs.
4. **Tuning-constant utilities:** `tuning_for_efficiency(family, eff)` and `tuning_for_breakdown(family, bdp)` — replicate `lmrob.control`'s logic for setting `tuning.psi` (95% efficiency) and `tuning.chi` (50% breakdown) from named arguments.

**Acceptance per family:**

```python
# pseudocode for tests/unit/test_psi.py
@pytest.mark.parametrize("family", ALL_FAMILIES)
@pytest.mark.parametrize("kind", ["rho", "psi", "psi_prime", "wgt"])
def test_psi_matches_r(r_session, family, kind):
    x = np.linspace(-5, 5, 1001)
    k = default_tuning(family)
    py = getattr(pylmrob.psi, kind)(x, family, k)
    r  = r_session.Mpsi(x, k, family) if kind == "psi" else ...
    assert_allclose(py, r, rtol=1e-12, atol=1e-14)
```

- All families pass at `rtol=1e-12`.
- `Epsi2`, `EDpsi` match R's `.Mwgt.psi1(...)$Epsi2` etc. at `rtol=1e-10`.
- Tuning constant tables match R bit-for-bit (CSV diff is empty).

**Validation:**
```bash
pytest tests/unit/test_psi.py -v
python -m pylmrob.benchmarks.psi  # ensure ≤ 2× R speed for vectorized eval
```

---

### Phase 3 — Robust M-Scale (3 days)

**Inputs:** Phase 2 (chi function).

**Algorithm:** Iteratively solve `mean(chi(r/σ)) = b0` for σ, where `b0` is the consistency constant (typically 0.5 for 50% bdp).

**Tasks:**

1. **`_core/_scale.pyx`:**
   - `m_scale(r[:], chi_family, k, b0, max_iter, tol, init_scale=NaN)` → `(scale, converged, iter_count)`
   - Initial estimate: MAD if `init_scale` is NaN.
   - Iteration: σ_{k+1} = σ_k · sqrt(mean(chi(r/σ_k, k)) / b0). (R uses this multiplicative form.)
   - Convergence: `|σ_{k+1} - σ_k| < tol · σ_k`.
2. **`scale.py`** — public `m_scale(r, family="bisquare", k=..., b0=0.5, ...)`.
3. **Edge cases:** all-zero residuals → return 0; near-zero scale → R's `zero.tol` semantics; non-finite residuals → raise.

**Acceptance:**
- Match R's `lmrob.mscale(r, control)` to `rtol=1e-9` on 50 random residual vectors.
- Iteration count within ±2 of R.

**Validation:** `pytest tests/unit/test_scale.py -v`.

---

### Phase 4 — Fast-S Resampling (THE BIG ONE — 2 weeks)

This is the speed-critical heart. Read Salibian-Barrera & Yohai (2006) and `src/lmrob.c::fast_s` carefully before starting.

**Algorithm sketch:**

```
input: X (n×p), y (n), control
1. for t = 1..nResample:
     a. draw a p-subset I_t (with non-singular check, retry up to mts times)
     b. β_t ← solve X[I_t] β = y[I_t]
     c. for k = 1..K (default 2):                # I-step refinement
          σ_t ← M-scale of residuals(β_t)
          β_t ← one IRWLS step using ψ at current σ_t
     d. keep β_t in best-of-`best_r` heap by σ_t
2. for each of the best_r candidates:
     refine to convergence (full IRWLS until rel_tol)
3. return the candidate with smallest σ
```

**Tasks (split per session):**

#### 4.1 Subsampling primitives
`_core/_fast_s.pyx`:
- `singular_safe_subsample(rng_state, n, p, mts) → indices[p]`: draws p indices, attempts QR, retries on rank-deficient. Match R's `mts` parameter (max-tries-singular, default 1000 in KS2014).
- Use `cython_lapack.dgeqrf` + condition-number check.

#### 4.2 Subset solver
- `solve_subset(X, y, indices, β_out)`: given a p-subset, solves the p×p system via `dgesv` (or QR via `dgels` if better-conditioned).
- Returns success flag.

#### 4.3 IRWLS step
- `irwls_step(X, y, β, σ, ψ_family, k, β_out)`: one weighted-least-squares update using weights `w_i = ψ(r_i/σ) / (r_i/σ)`. Use `dgels` on `diag(√w) X` and `diag(√w) y`. Match R's handling of zero weights.

#### 4.4 Candidate refinement
- `refine_candidate(X, y, β, σ, control, β_out, σ_out, conv_out)`: full IRWLS to convergence.

#### 4.5 Best-r heap
- Fixed-size sorted insertion (best_r ≤ 20 typically). No need for a real heap; linear scan is fine.

#### 4.6 Outer driver
- `fast_s(X, y, control, seed, β_out, σ_out, info_out)`:
  - Allocates workspaces once.
  - `prange(nResample, nogil=True)` over the candidate loop.
  - Each thread has its own RNG (PCG64 seeded from a `SeedSequence` per thread).
  - Aggregates best_r results across threads (critical section).
  - Refines best_r in parallel.
  - Returns best.

#### 4.7 Large-n branch
- Match R's `fast_s_large_n`: when `n > eps.x * p`, partition data into `groups` random groups of size `n_group`, run fast-S on each, take best β per group, then refine across full data. Match defaults from R (`groups=5`, `n.group=400`, threshold `1.5 * p`).

#### 4.8 Python wrapper
- `_core/__init__.py` exposes `fast_s(X, y, control)` returning `(beta, scale, info)` named tuple.

**Acceptance:**

- For 30+ datasets and 10 different seeds each: scale `σ` matches R's `init.S$scale` to `rtol=1e-5` (occasional outliers acceptable due to RNG, document each).
- Coefficients match to `rtol=1e-4` post-S (will be tightened after MM in Phase 6).
- For `n=10000, p=20`, single-threaded: ≤ 1.3× R wall-clock. Multi-threaded with 4 cores: ≤ 0.5× R.
- `cython -a` HTML inspection: zero yellow lines inside the resampling loop.

**Validation:**
```bash
pytest tests/unit/test_fast_s.py -v
python tests/benchmarks/bench_fast_s.py --compare-r
```

---

### Phase 5 — M-S Estimator for Factor Designs (1 week)

**Why:** Pure-S resampling fails on designs with categorical predictors because random p-subsets often produce singular subdesigns. Maronna & Yohai (2000) propose M-S: separate continuous from categorical part, alternate.

**Tasks:**

1. Implement `split_design(X, factor_cols)` in `_utils.py`: separates `X = [X_cat | X_cont]`.
2. **`ms_estimator.py`:**
   - `m_s_fit(X_cat, X_cont, y, control)`:
     - Initialize β_cat by L1 fit on factor part.
     - Run S on continuous part of partial residuals.
     - Iterate: update β_cat by weighted L1, β_cont by S/MM.
     - Use `scipy.optimize.linprog` or a direct simplex for the L1 sub-problem (validate against R's `quantreg::rq` if needed).
3. Auto-detection in `lmrob.py`: if any column of X is categorical (formulaic provides this metadata), default `init="M-S"`.

**Acceptance:**
- Match R on `education`, `pension`, and 3+ synthetic factor-heavy designs.
- `init="S"` and `init="M-S"` paths both produce results matching R.

**Validation:** `pytest tests/integration/test_factor_designs.py -v`.

---

### Phase 6 — MM Iteration (3–4 days)

**Inputs:** Phase 4/5 (initial S or M-S estimate).

**Algorithm:** Given (β_S, σ_S), iterate IRWLS using the higher-efficiency ψ (typically tuned for 95% efficiency) holding σ_S fixed, until β converges.

**Tasks:**

1. **`_core/_mm.pyx::mm_iterate(X, y, β_init, σ, control)`** — IRWLS loop with the efficiency-tuned ψ. Convergence on `‖β_{k+1} - β_k‖ / ‖β_k‖ < rel_tol`.
2. **Step-halving** when residual sum-of-squares of weighted residuals increases — R does this; match it.
3. **Iteration cap** with `converged` flag in result.

**Acceptance:**
- Final β matches R to `rtol=1e-6` on the full validation corpus.
- Final σ matches to `rtol=1e-8` (σ doesn't change in MM).
- Final weights match to `rtol=1e-6`.

**Validation:** `pytest tests/integration/test_lmrob_api.py::test_mm_matches_r -v`.

---

### Phase 7 — Inference & Covariance (4 days)

**Three estimators in R:**

- `.vcov.avar1` (default for KS2014): asymptotic variance from sandwich form.
- `.vcov.w` (Koller-Stahel 2011): adjusted for small samples.
- `Asymp` (legacy MM): older form.

**Tasks:**

1. `inference.py::vcov_avar1(X, weights, ψ, σ, ...)`: implement sandwich `(X'WX)^-1 X'W²X (X'WX)^-1` with the appropriate scaling factors `1/(n-p) · Σψ²(r/σ) / (Σψ'(r/σ))² · σ²`.
2. `inference.py::vcov_w(...)`: KS2011 finite-sample corrections.
3. `inference.py::vcov_asymp(...)`: legacy.
4. `LmRobResults.summary()` — t-statistics, p-values (uses `scipy.stats.t` or `norm` per R's choice), CIs.
5. Wald test, `predict(new_X, interval="confidence"|"prediction")`.

**Acceptance:**
- `cov` matches R to `rtol=1e-6` on the corpus.
- t-stats and p-values match.
- `summary()` output reproduces R's column structure.

**Validation:** `pytest tests/unit/test_inference.py -v`.

---

### Phase 8 — Public API & Control (3 days)

**Tasks:**

1. **`control.py::Control` dataclass** with every R `lmrob.control` parameter. Provide `Control.preset("KS2014")`, `Control.preset("KS2011")`, `Control.preset("MM")`. Defaults must match R **exactly**, by setting.
2. **`lmrob(formula, data, control=None, weights=None, na_action="drop", ...)`** — top-level fit function. Returns `LmRobResults`.
3. **scikit-learn-style class:** `LmRob(control=...).fit(X, y).predict(X_new)`.
4. **`formula.py`** — `formulaic` integration; handle `C(x)`, `I(x**2)`, factors, `0+`, etc. Test against R's `model.matrix` output.
5. **`results.py::LmRobResults`** — properties: `coef_`, `scale_`, `weights_`, `residuals_`, `fitted_`, `cov_`, `df_residual_`, `converged_`, `init_`, `rweights_`, `psi_`, `nobs_`. Methods: `summary()`, `predict()`, `confint()`, `__repr__()`.

**Acceptance:**

- `lmrob("y ~ x1 + x2", df)` works end-to-end.
- `repr()` and `summary()` resemble R's output text.
- Round-trip pickle works.

**Validation:** `pytest tests/integration/test_lmrob_api.py -v`.

---

### Phase 9 — Diagnostics (3 days)

**Tasks:**

1. `diagnostics.py::plot(results)` — 4-panel: (a) Residuals vs Fitted, (b) Normal Q-Q of standardized residuals, (c) Robust distance vs Residual (the "Robust diagnostic plot"), (d) Residuals vs Leverage. Match R's `plot.lmrob`.
2. `cooks_distance(results, robust=True)`.
3. Influence statistics: `hatvalues` (robust), `dfbetas`-like.

**Acceptance:** Visual snapshot tests with `pytest-mpl`. Diagnostic statistics match R `lmrob`'s `weights`, `rweights` to `rtol=1e-6`.

---

### Phase 10 — Validation & Property Tests (parallel with all phases; 1 week to consolidate)

**Tasks:**

1. **Classical-data sweep** (`test_vs_r_classical.py`): all 11 datasets × 4 settings → JSON snapshots assert.
2. **Synthetic sweep** (`test_vs_r_synthetic.py`): n × p × contamination grid; assert convergence rate ≥ R's, β within tolerance.
3. **Property tests** (`test_equivariance.py` with Hypothesis):
   - **Affine equivariance:** `lmrob(y, X)` then `lmrob(a + b·y, X·M)` should give predictably-transformed coefficients.
   - **Regression equivariance:** adding linear function of X to y shifts β by that linear function.
   - **Scale equivariance:** scaling X scales β inversely.
4. **Edge-case suite:** all-equal y, perfect fit, exactly p observations, n < p (must error gracefully), NaN/Inf handling.
5. **RNG reproducibility:** same seed → same output across runs and platforms (within `rtol=1e-12`).

**Acceptance:**

- 100% of corpus tests pass at agreed tolerances.
- Hypothesis runs 500 examples per property test in CI without finding violations.

---

### Phase 11 — Performance (1 week)

**Tasks:**

1. **Benchmark harness** `tests/benchmarks/bench_full_fit.py`: matrix of (n, p, settings) measured against R using `microbenchmark` via rpy2. Output as `bench-report.md` committed per release.
2. **Profile** with `py-spy` and `cython -a`. Targets:
   - Inner subsample loop: zero Python overhead.
   - LAPACK calls: direct `cython_lapack`, no SciPy Python wrappers.
   - IRWLS: BLAS-3 via `dsyrk` for `X'WX`.
3. **Memory:** preallocate workspaces; use `np.empty` not `np.zeros`; avoid temporary arrays in hot path.
4. **Parallelism tuning:** OpenMP chunk size, dynamic vs static scheduling. Match R's behavior under `OMP_NUM_THREADS`.
5. **Numerical stability vs speed tradeoffs:** document where we choose `dgels` (QR) over `dposv` (Cholesky) for stability.

**Acceptance (single-threaded, on a modern x86_64 with OpenBLAS):**

| n | p | R time | Python target |
|---|---|---|---|
| 100 | 5 | 50 ms | ≤ 75 ms |
| 1,000 | 10 | 200 ms | ≤ 250 ms |
| 10,000 | 20 | 5 s | ≤ 5 s |
| 100,000 | 50 | 90 s | ≤ 90 s |

Multi-threaded (4 cores): each row should be 2–3× faster; R does not parallelize as cleanly.

---

### Phase 12 — Packaging, Docs, Release (4 days)

**Tasks:**

1. **Sphinx docs** at `docs/`: API reference (autodoc), tutorial mirroring `lmrob` vignette, "Migrating from R" guide, "Comparison with statsmodels.RLM" section.
2. **README** with: install, quickstart, comparison table to alternatives, citation, license.
3. **`cibuildwheel` matrix:** `cp310-*`, `cp311-*`, `cp312-*`, `cp313-*` for `manylinux_2_28_x86_64`, `manylinux_2_28_aarch64`, `macosx_x86_64`, `macosx_arm64`, `win_amd64`. ARM64 cross-build via QEMU on Linux.
4. **PyPI release** via trusted publisher (OIDC).
5. **Conda-forge feedstock** PR (post-PyPI).
6. **Zenodo DOI** for citation.

**Acceptance:**

- `pip install pylmrob` works on all target platforms in clean venvs.
- `import pylmrob; pylmrob.lmrob(...)` works without compilation on user machine.
- Docs build clean.

---

## 5. Cross-cutting Concerns

### 5.1 Numerical Tolerances

Default tolerances in `tests/conftest.py`:

| Output | rtol | atol | Notes |
|---|---|---|---|
| psi/chi/wgt | 1e-12 | 1e-14 | Pure numeric |
| M-scale | 1e-9 | 1e-12 | Iterative |
| Init S β | 1e-4 | 1e-6 | RNG-dependent |
| Init S σ | 1e-5 | 1e-8 | RNG-dependent |
| MM β | 1e-6 | 1e-8 | After convergence |
| MM σ | 1e-8 | 1e-10 | Inherited from S |
| Weights | 1e-6 | 1e-8 | |
| Cov | 1e-6 | 1e-8 | |

Tolerances may be loosened **only with a comment in the test explaining why**.

### 5.2 RNG Strategy

- Seed accepted as `int | np.random.Generator | None` in `lmrob(seed=...)`.
- Internal: `SeedSequence(seed).spawn(num_threads)` → one PCG64 per thread.
- Document explicitly: "Results are bit-reproducible across runs with the same seed and thread count, but **not** identical to R's `set.seed(...)` results because R uses Mersenne Twister."

### 5.3 Error Handling

- Singular full design → `LinAlgError` with a clear message naming the rank-deficient columns.
- Convergence failure → return result with `converged=False` and a `RuntimeWarning`. Do **not** raise; R doesn't.
- All-zero scale (perfect fit) → match R's behavior (return σ=0, weights=1, warn).

### 5.4 Documentation Requirements per Module

Every public function: NumPy-style docstring with Examples block that runs under `pytest --doctest-modules`. Every C/Cython function: Doxygen-style comment block referencing the source paper section.

### 5.5 What to do when results disagree with R

Process:
1. Reproduce in minimal form. Save inputs + seeds.
2. Run R with `trace.lev = 4`; run Python with the equivalent verbose mode (add one if needed).
3. Diff iteration-by-iteration logs.
4. Either fix the bug or document the divergence in `docs/numerical-notes.md` with rationale.

---

## 6. Risks & Mitigations

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Edge cases in fast-S (singular subsamples, scale collapse) | High | High | Port R's accumulated guards verbatim into Cython; large randomized test corpus |
| `ggw`/`lqq` tuning constants wrong | Medium | High | Generate via R script and embed; validate by exact-equality test |
| OpenMP issues on macOS | Medium | Medium | Document `brew install libomp`; provide single-threaded fallback build |
| Cython ABI breakage with new NumPy | Low | Medium | Pin via `oldest-supported-numpy` in build deps |
| RNG divergence makes validation noisy | High | Low | Use `nResample` large enough that S-step finds same basin >99% of time; mark RNG-sensitive tests with `@pytest.mark.rng_sensitive` and run k seeds |
| GPL licensing limits adoption | Medium | Medium | Document clearly; recommend `statsmodels.RLM` for users who can't use GPL |
| `formulaic` API changes | Low | Low | Pin minor version; abstract behind `formula.py` |
| Sole-developer bus factor | High | High | Detailed PLAN.md (this file), thorough docs, recorded design decisions |

---

## 7. Definition of Done (v1.0)

- [ ] All 12 phases complete.
- [ ] All tests pass on Linux/macOS/Windows × Python 3.10–3.13.
- [ ] Validation corpus matches R within documented tolerances.
- [ ] Performance targets met on the reference machine (record machine specs in `bench-report.md`).
- [ ] Wheels published to PyPI for all target platforms.
- [ ] Docs published to GitHub Pages.
- [ ] README has install + quickstart + citation.
- [ ] CHANGELOG up to date.
- [ ] Tagged `v1.0.0` with release notes.
- [ ] Conda-forge feedstock submitted.
- [ ] At least one external user has built against it (sanity check).

---

## 8. Working Style for Claude Code

- **Branch per task.** One PR per task. Keep PRs ≤ 500 lines of diff where possible.
- **Tests before / with code.** Every numeric task has its R-reference test landing in the same PR as the implementation.
- **Pre-commit.** Ruff format + lint, mypy strict, cython linting must all pass.
- **Commit messages:** Conventional Commits (`feat:`, `fix:`, `perf:`, `test:`, `docs:`, `build:`).
- **When stuck on numerical disagreement:** dump iteration logs from both R and Python, diff line by line. Don't guess.
- **When tempted to hand-write LAPACK:** stop, use `cython_lapack`.
- **When tempted to skip a test:** stop, write the test.
- **Document every deviation from R in `docs/numerical-notes.md`** with a code reference and rationale.

---

## 9. Quick Reference: Key R Symbols → Python Targets

| R | Python | Phase |
|---|---|---|
| `lmrob(formula, data, ...)` | `pylmrob.lmrob(formula, data, ...)` | 8 |
| `lmrob.control(...)` | `pylmrob.Control(...)` | 8 |
| `lmrob.S(x, y, control)` | `pylmrob._core.fast_s(X, y, control)` | 4 |
| `lmrob..M..fit(...)` | `pylmrob._core.mm_iterate(...)` | 6 |
| `lmrob.M.S(x, y, control)` | `pylmrob.ms_estimator.m_s_fit(...)` | 5 |
| `Mpsi(x, c, "bisquare")` | `pylmrob.psi.psi(x, "bisquare", c)` | 2 |
| `Mchi(x, c, "bisquare")` | `pylmrob.psi.rho(x, "bisquare", c)` | 2 |
| `Mwgt(x, c, "bisquare")` | `pylmrob.psi.wgt(x, "bisquare", c)` | 2 |
| `lmrob.mscale(r, control)` | `pylmrob.scale.m_scale(r, ...)` | 3 |
| `.vcov.avar1`, `.vcov.w` | `pylmrob.inference.vcov_*` | 7 |
| `summary.lmrob(fit)` | `fit.summary()` | 8 |
| `plot.lmrob(fit)` | `pylmrob.diagnostics.plot(fit)` | 9 |

---

## 10. Estimated Timeline

| Phase | Solo dev (calendar weeks) |
|---|---|
| 0. Bootstrap | 0.3 |
| 1. R harness | 0.5 |
| 2. Psi functions | 1.0 |
| 3. M-scale | 0.5 |
| 4. Fast-S | 2.0 |
| 5. M-S init | 1.0 |
| 6. MM iteration | 0.6 |
| 7. Inference | 0.8 |
| 8. Public API | 0.6 |
| 9. Diagnostics | 0.6 |
| 10. Validation consolidation | 1.0 |
| 11. Performance | 1.0 |
| 12. Packaging & release | 0.6 |
| **Total** | **~10.5 weeks** |

Add ~30% buffer for unknown unknowns → ~14 weeks to v1.0.

---

*End of PLAN.md.*
