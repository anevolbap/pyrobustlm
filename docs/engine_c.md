# The `engine_c=True` fast path

`Control(engine_c=True)` is the default since v0.5.11. The fit runs
through a monolithic Cython kernel that does the entire S-init,
survivor refinement, MM iteration, optional D-scale
(Koller & Stahel 2014), and `vcov_avar1` inside one `nogil` C block
with a single workspace allocation.

The kernel mirrors the structure of `robustbase/src/lmrob.c`:
- LAPACK via `scipy.linalg.cython_lapack` (`dgesv` for the p-subset
  solve and the `(X'WX)` inversion, `dgels` for IRWLS, `dsyev` for the
  vcov posdefify).
- Subset draws via numpy's `bitgen_t` capsule (Floyd's combination
  algorithm).
- Per-family `psi` / `psi'` / `chi` / `chi'` kernels inlined as `cdef`
  helpers; family dispatch is a small constant-folded enum.

## When to use it

Turn it on when you want the fastest end-to-end fit and you can
tolerate basin-of-attraction drift (see below):

```python
from pylmrob import lmrob, Control

fit = lmrob(formula, df, control=Control(engine_c=True), seed=42)
```

It currently covers all five lmrob-supported psi families
(`bisquare`, `hampel`, `optimal`, `lqq`, `ggw`) and both
`cov=".vcov.avar1"` and the `setting="KS2014"` / `"KS2011"` SMDM
pipeline.

## Behaviour at larger n

`engine_c=True` runs the monolithic Cython kernel only when
`n*p^2 < 100,000`. Above that, the threaded default path is faster
because the monolithic kernel is a single non-parallel call, whereas
the default `_resample_chunk` parallelises across a thread pool.

The fallback is automatic: setting `Control(engine_c=True)` also
enables `n_workers=0` (auto threading) on the default path when the
problem is large. Users do not need to choose a path manually.

Concrete example at n=5000, p=30, nResample=500 (single-thread
OpenBLAS):

| config | wall-clock |
|---|---|
| `Control()` (default, serial) | ~2200 ms |
| `Control(n_workers=0)` | ~860 ms |
| `Control(engine_c=True)` | ~870 ms (auto fallback) |

## What you give up

The Cython subset-draw uses Floyd's combination algorithm via numpy's
`bitgen_t`. That's mathematically equivalent to
`np.random.Generator.choice(n, p, replace=False)` but the byte sequence
consumed from the BitGenerator differs, so the basin of attraction can
shift slightly. On tiny-n problems (e.g. n=21 stackloss with some
seeds) the shifted basin occasionally lands on a fit where `X'WX` is
singular and `vcov_avar1` fails. `lmrob()` catches that
`FloatingPointError` once and retries with `engine_c=False` so the
fit always succeeds. The fit's coef/scale stay within ~1e-3 rerr of
the numpy path on those edge cases. To get byte-identical fits with
the pre-engine_c releases, pass `Control(engine_c=False)`.

Note that `Control(rng="R")` forces `engine_c=False` automatically: the
monolithic engine owns its own BitGenerator and can't share R-mode's
``unif_rand`` state.

## Benchmark

Single-threaded BLAS, 7-rep median, four classical datasets × three
settings, on a 16-core OpenBLAS Linux x86_64 box. R wall-clock on the
same fits is 3-7 ms.

| dataset / setting | default | engine_c | speedup | vs R |
|---|---|---|---|---|
| stackloss MM     | 22.4 ms |  4.0 ms |  5.6x | 1.3x R |
| stackloss KS2014 | 62.8 ms |  7.1 ms |  8.8x | 1.4x R |
| stackloss KS2011 | 70.4 ms |  7.2 ms |  9.8x | 1.4x R |
| delivery MM      | 24.7 ms |  5.1 ms |  4.8x | 1.7x R |
| delivery KS2014  | 68.8 ms |  8.7 ms |  7.9x | 1.7x R |
| phosphor MM      | 23.5 ms |  3.7 ms |  6.4x | 1.2x R |
| phosphor KS2014  | 69.5 ms |  6.1 ms | 11.4x | 1.2x R |
| phosphor KS2011  | 74.2 ms |  6.2 ms | 12.0x | 1.2x R |
| salinity MM      | 26.8 ms |  5.2 ms |  5.2x | 1.7x R |
| salinity KS2014  | 76.0 ms |  9.9 ms |  7.7x | 1.7x R |

Reproduce with `scripts/bench_engine_c.py`.

## What's still pure Python

The remaining Python work per fit (with the default `engine_c=True`) is:
- the formula parser (about 0.3 ms via a hand-written fast path for
  the common `y ~ x1 + x2 + ...` pattern; complex formulas go through
  `formulaic`)
- `Control` field unpacking and a few small ndarray allocations
- the final `LmRobResults` construction

The user-facing array API (`LmRob.fit(X, y)`) skips formula parsing
entirely. The total wall-clock floor is approximately the underlying
BLAS work (the same R does internally) plus one Python/C boundary
cross. We pay that crossing because `pylmrob` keeps a Python
public API; R can stay in C from `R_lmrob_S` all the way to the
returned `SEXP`.
