# R source map

Mapping `robustbase`'s R and C source to the corresponding Python module
in `pylmrob`. Update this as each phase lands.

## R sources (CRAN mirror `cran/robustbase`, package `R/`)

| R file | What it is | Python target | Phase |
|--------|------------|---------------|-------|
| `R/lmrob.R` | Top-level `lmrob()` and dispatch | `pylmrob.lmrob` | 8 |
| `R/lmrob.MM.R` | MM iteration | `pylmrob._core._mm` + `pylmrob.lmrob` | 6 |
| `R/lmrob.M.S.R` | M-S initial | `pylmrob.ms_estimator` | 5 |
| `R/lmrobMMcontrol.R` | `lmrob.control` defaults; psi-tuning constant tables | `pylmrob.control` + `pylmrob._core._ctables` | 8, 2 |
| `R/lmrobMMresid.R` | residuals, predict, deviance | `pylmrob.results` | 8 |
| `R/Mpsi.R`, `R/Mchi.R`, `R/Mwgt.R` | psi/chi/wgt API and tuning helpers | `pylmrob.psi` (wraps `pylmrob._core._psi`) | 2 |
| `R/lmrob.mscale.R` | M-scale wrapper | `pylmrob.scale` (wraps `pylmrob._core._scale`) | 3 |
| `R/.vcov.avar1.R`, `R/.vcov.w.R` | covariance estimators | `pylmrob.inference` | 7 |

## C sources (`src/`)

| C file | What it is | Python target | Phase |
|--------|------------|---------------|-------|
| `src/lmrob.h` | constants, structs | reference; types live in `_core/_*.pxd` if needed | various |
| `src/lmrob.c::fast_s` | core resampling loop | `pylmrob._core._fast_s::fast_s` | 4 |
| `src/lmrob.c::fast_s_large_n` | large-n branch | `pylmrob._core._fast_s::fast_s_large_n` | 4.7 |
| `src/lmrob.c::rwls` | IRWLS step | `pylmrob._core._fast_s::irwls_step` and `pylmrob._core._mm::mm_iterate` | 4, 6 |
| `src/lmrob.c::rho` etc. | psi/chi inner kernels | `pylmrob._core._psi` | 2 |
| `src/lmrob-psifuns.c` | psi/chi/wgt/Epsi for all 6 families | `pylmrob._core._psi` | 2 |
| `src/mc.c` | Mahalanobis-like utilities | reference only | n/a |
| `src/init.c` | `.Call` registrations | reference only | n/a |

## Symbols cross-reference

| R symbol | Python target | Phase |
|----------|---------------|-------|
| `lmrob(formula, data, ...)` | `pylmrob.lmrob(formula, data, ...)` | 8 |
| `lmrob.control(...)` | `pylmrob.Control(...)` / `Control.preset(...)` | 8 |
| `lmrob.S(x, y, control)` | `pylmrob._core._fast_s.fast_s(X, y, control)` | 4 |
| `lmrob..M..fit(...)` | `pylmrob._core._mm.mm_iterate(...)` | 6 |
| `lmrob.M.S(x, y, control)` | `pylmrob.ms_estimator.m_s_fit(...)` | 5 |
| `Mpsi(x, c, "bisquare")` | `pylmrob.psi.psi(x, "bisquare", c)` | 2 |
| `Mchi(x, c, "bisquare")` | `pylmrob.psi.rho(x, "bisquare", c)` | 2 |
| `Mwgt(x, c, "bisquare")` | `pylmrob.psi.wgt(x, "bisquare", c)` | 2 |
| `lmrob.mscale(r, control)` | `pylmrob.scale.m_scale(r, ...)` | 3 |
| `.vcov.avar1`, `.vcov.w` | `pylmrob.inference.vcov_avar1`, `vcov_w` | 7 |
| `summary.lmrob(fit)` | `LmRobResults.summary()` | 8 |
| `plot.lmrob(fit)` | `pylmrob.diagnostics.plot(fit)` | 9 |
