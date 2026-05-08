# R source map

Mapping `robustbase`'s R and C source to the corresponding Python module
in `pyrobustlm`. Update this as each phase lands.

## R sources (CRAN mirror `cran/robustbase`, package `R/`)

| R file | What it is | Python target | Phase |
|--------|------------|---------------|-------|
| `R/lmrob.R` | Top-level `lmrob()` and dispatch | `pyrobustlm.lmrob` | 8 |
| `R/lmrob.MM.R` | MM iteration | `pyrobustlm._core._mm` + `pyrobustlm.lmrob` | 6 |
| `R/lmrob.M.S.R` | M-S initial | `pyrobustlm.ms_estimator` | 5 |
| `R/lmrobMMcontrol.R` | `lmrob.control` defaults; psi-tuning constant tables | `pyrobustlm.control` + `pyrobustlm._core._ctables` | 8, 2 |
| `R/lmrobMMresid.R` | residuals, predict, deviance | `pyrobustlm.results` | 8 |
| `R/Mpsi.R`, `R/Mchi.R`, `R/Mwgt.R` | psi/chi/wgt API and tuning helpers | `pyrobustlm.psi` (wraps `pyrobustlm._core._psi`) | 2 |
| `R/lmrob.mscale.R` | M-scale wrapper | `pyrobustlm.scale` (wraps `pyrobustlm._core._scale`) | 3 |
| `R/.vcov.avar1.R`, `R/.vcov.w.R` | covariance estimators | `pyrobustlm.inference` | 7 |

## C sources (`src/`)

| C file | What it is | Python target | Phase |
|--------|------------|---------------|-------|
| `src/lmrob.h` | constants, structs | reference; types live in `_core/_*.pxd` if needed | various |
| `src/lmrob.c::fast_s` | core resampling loop | `pyrobustlm._core._fast_s::fast_s` | 4 |
| `src/lmrob.c::fast_s_large_n` | large-n branch | `pyrobustlm._core._fast_s::fast_s_large_n` | 4.7 |
| `src/lmrob.c::rwls` | IRWLS step | `pyrobustlm._core._fast_s::irwls_step` and `pyrobustlm._core._mm::mm_iterate` | 4, 6 |
| `src/lmrob.c::rho` etc. | psi/chi inner kernels | `pyrobustlm._core._psi` | 2 |
| `src/lmrob-psifuns.c` | psi/chi/wgt/Epsi for all 6 families | `pyrobustlm._core._psi` | 2 |
| `src/mc.c` | Mahalanobis-like utilities | reference only | n/a |
| `src/init.c` | `.Call` registrations | reference only | n/a |

## Symbols cross-reference

| R symbol | Python target | Phase |
|----------|---------------|-------|
| `lmrob(formula, data, ...)` | `pyrobustlm.lmrob(formula, data, ...)` | 8 |
| `lmrob.control(...)` | `pyrobustlm.Control(...)` / `Control.preset(...)` | 8 |
| `lmrob.S(x, y, control)` | `pyrobustlm._core._fast_s.fast_s(X, y, control)` | 4 |
| `lmrob..M..fit(...)` | `pyrobustlm._core._mm.mm_iterate(...)` | 6 |
| `lmrob.M.S(x, y, control)` | `pyrobustlm.ms_estimator.m_s_fit(...)` | 5 |
| `Mpsi(x, c, "bisquare")` | `pyrobustlm.psi.psi(x, "bisquare", c)` | 2 |
| `Mchi(x, c, "bisquare")` | `pyrobustlm.psi.rho(x, "bisquare", c)` | 2 |
| `Mwgt(x, c, "bisquare")` | `pyrobustlm.psi.wgt(x, "bisquare", c)` | 2 |
| `lmrob.mscale(r, control)` | `pyrobustlm.scale.m_scale(r, ...)` | 3 |
| `.vcov.avar1`, `.vcov.w` | `pyrobustlm.inference.vcov_avar1`, `vcov_w` | 7 |
| `summary.lmrob(fit)` | `LmRobResults.summary()` | 8 |
| `plot.lmrob(fit)` | `pyrobustlm.diagnostics.plot(fit)` | 9 |
