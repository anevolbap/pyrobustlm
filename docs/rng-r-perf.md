# `Control(rng="R")` performance

`Control(rng="R")` drives the fast-S resample loop through `pylmrob.r_set_seed` +
`r_sample_noreplace` / `r_subsample_nonsingular` so subset draws are byte-identical
to robustbase. This page benchmarks the wall-clock cost vs the default PCG64 path.

Companion of [`bench-report.md`](bench-report.md), which compares pylmrob (PCG64,
default) against R itself.

## Setup

- pylmrob 0.5.16
- Cython R RNG kernel enabled (`pylmrob._core._r_rng`)
- Python 3.13.12, NumPy + scipy with system OpenBLAS, Linux x86_64
- 5 reps per case, median reported

## Wall-clock: PCG64 vs R-mode

| case          | n    | p  | PCG64 (ms) | R-mode (ms) | R-mode / PCG64 |
|---------------|-----:|---:|-----------:|------------:|---------------:|
| stackloss     |   21 |  4 |       27.3 |         9.7 |          0.36x |
| synth_500_5   |  500 |  5 |       42.9 |        80.1 |          1.87x |
| synth_2000_10 | 2000 | 10 |      246.7 |       389.9 |          1.58x |

`Control(rng="R")` is faster than PCG64 on small problems (the formula parser and
PCG64 BitGenerator spinup amortise better at scale) and slower on larger problems
by ~1.6-1.9x. The slowdown comes from:

- Forced `n_workers=1` (R-mode is sequential; no thread-pool benefit).
- Forced `engine_c=False` (the monolithic Cython engine owns its BitGenerator and
  can't share R-mode state).
- The R-style subset draw is a full permutation per attempt (n `unif_rand` draws);
  the default Cython Floyd's combination does `p` per attempt.

## Microbench: draw helpers

The Cython R-RNG kernel (`pylmrob._core._r_rng`) gives substantial speedups over
the pure-Python fallback for the standalone helpers:

| op                                                  | pure Python | Cython | speedup |
|-----------------------------------------------------|------------:|-------:|--------:|
| `r_sample_noreplace(n=200, k=200)` × 500            |     133 ms | 1.0 ms |   ~130x |
| `r_subsample_nonsingular(n=200, p=5)` × 500         |     162 ms | 1.5 ms |   ~108x |

## When to use `rng="R"`

- Reproducing a fit from an R analysis where someone published `set.seed(seed)`.
- Generating reference data for cross-language validation.
- Tighter agreement with R's `lmrob` (rtol ~1.7e-5 on coefs vs ~1e-3 for PCG64
  on the basin-drift cases — see [`numerical-notes.md`](numerical-notes.md)).

Not recommended for hot loops on large data: the default `rng="PCG64"` with
threading is faster and statistically just as valid.

## Remaining work

The residual rtol~1e-5 gap to R's `lmrob` lives in the resample loop's
"associated scale" computation (R uses `find_scale`; pylmrob uses
`_mscale_generic` in the Cython kernel). Porting `refine_fast_s` + `find_scale`
to Cython end-to-end is tracked in [`../plan.md`](../plan.md#11-post-v0516-r-compat-track).
