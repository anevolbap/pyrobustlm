# Parallelising the `engine_c` kernel

Design notes for splitting `cy_lmrob_fit` / `cy_lmrob_fast_s` across
multiple threads. Today the monolithic engine is one `nogil` block
that does not parallelise internally; `lmrob()` auto-falls-back to the
threaded NumPy path when `n*p^2 >= 100k`. If the kernel itself could
use multiple threads, the fallback heuristic could go away and the
median bench `engine_c/R` ratio (currently 0.93x, with the worst
large-n case at 0.62x R) could drop further.

## Where the parallel opportunity lives

In `src/pyrobustlm/_core/_lmrob.pyx`, inside `cy_lmrob_fast_s`:

    for try_i in range(nResample):
        # draw a p-subset (Floyd via bitgen)
        # initial solve (dgesv on the p-subset)
        # k_fast_s K-step refinements (one m-scale + IRWLS each)
        # maintain best-of-best_r heap

Each iteration is independent of the others except for the shared
`best_*` heap. The work is dominated by the K-step refinements (m-scale
+ dgels-based IRWLS at full n). At `nResample=500` and `n=5000, p=30`
the loop is ~700 ms; trivially embarrassingly parallel.

## Constraints

1. **Per-thread scratch.** `_Scratch` holds `sub_X`, `sub_y`, `X_w`,
   `y_w`, `r`, `w`, `perm`, `beta`, `beta_prev`, LAPACK `work` and
   `ipiv`. Today it's one struct allocated at the top of
   `cy_lmrob_fast_s`. For N threads we need N copies. Allocation cost
   is O(N * n * p) bytes (~2 MB at n=5000 p=30 across 8 threads), which
   is fine.

2. **Bitgen sharing.** Numpy's `bitgen_t` C API is not thread-safe;
   concurrent `next_uint64` calls on the same bitgen will race. Three
   options:

   - **One bitgen per thread.** Spawn N child bitgens before the
     parallel region, each seeded with a thread-id-dependent stream
     from the user's seed (e.g. `np.random.SeedSequence(seed).spawn(N)`).
     Pre-allocate N bitgen capsules and pass an array of pointers into
     the kernel. Best determinism: same seed always gives the same
     fit independent of how many threads ran, as long as the splitting
     is deterministic.
   - **Lock around bitgen calls.** Cheaper to implement but throws
     away most of the parallel win (subset-draw is a hot inner loop
     called `mts` times per candidate).
   - **Per-thread reservoirs.** Pre-draw all needed uint64 values
     into a per-thread buffer before the parallel region. Same
     determinism as option 1 but more memory and more code.

   **Pick:** option 1. The spawn-from-SeedSequence pattern is what
   `numpy.random.Generator.spawn()` is for.

3. **Best-of-best_r heap.** Shared across all candidates. Two
   options:

   - **Per-thread heap + final merge.** Each thread keeps its own
     best-of-best_r/N heap, after the parallel region we merge the
     N heaps. Lock-free during the loop; merge is O(N * best_r) which
     is negligible.
   - **Critical-section update.** Lock around the heap update. Simpler
     code but a contention point on small problems.

   **Pick:** per-thread heap + merge. Same fit selection as the
   serial code (the merge sorts identically), no locking on the hot
   path.

4. **Determinism.** With option 1 + per-thread heap merge, the fit is
   deterministic given `(seed, n_workers)`. It is NOT deterministic
   across `n_workers` values, because the per-thread RNG streams
   differ. That is the same trade-off the existing `_fast_s` threaded
   path already accepts, and it matches the way `lmrob.engine_c` users
   already think about reproducibility (the comment in `Control` calls
   out RNG byte-level non-portability).

5. **OpenMP availability.** `meson.build` already has
   `omp_dep = dependency('openmp', required: false)` — best-effort,
   no-op on platforms without it. The kernel falls back to serial on
   non-OMP builds.

## API shape

`cy_lmrob_fit` and `cy_lmrob_fast_s` gain an `int n_workers` argument
(0 = auto, 1 = serial, >1 = explicit). `lmrob()` forwards
`control.n_workers` (which already exists for the NumPy path).

Inside the kernel:

    cdef int N = n_workers if n_workers > 0 else <heuristic>(n, p)
    cdef _Scratch* scratches = <_Scratch*>malloc(N * sizeof(_Scratch))
    for i in range(N):
        _alloc_scratch(&scratches[i], n, p)
    cdef bitgen_t** bgs = <bitgen_t**>malloc(N * sizeof(bitgen_t*))
    for i in range(N):
        bgs[i] = ...  # from the spawned capsules
    cdef double* best_scales = <double*>malloc(N * best_r * sizeof(double))
    cdef double* best_betas  = <double*>malloc(N * best_r * p * sizeof(double))

    with nogil, parallel(num_threads=N):
        tid = openmp.omp_get_thread_num()
        for try_i in prange(nResample, schedule="static"):
            # use scratches[tid], bgs[tid]
            # update best_scales[tid*best_r:..], best_betas[tid*best_r:..]

    # serial merge of N best heaps into the global best_r
    # survivor refinement on the global winner

## Smallest landable first step

Two PRs:

1. **Refactor `cy_lmrob_fast_s` to take an explicit scratch + bitgen
   pointer pair**, callable in a loop. Keep the existing public
   signature; introduce `_cy_lmrob_fast_s_chunk(start, end, scr,
   bg, best_scales, best_betas, ...)` as an internal helper. No
   parallelism yet. Tests pass unchanged.

2. **Add `n_workers` parameter** to the public kernel. Allocate N
   scratches and N bitgen capsules at the Python boundary (using
   `np.random.SeedSequence.spawn`). Use `cython.parallel.prange`
   over the chunk loop. Merge per-thread heaps. Update
   `lmrob.py` to forward `control.n_workers` and remove the
   `_engine_c_too_big` auto-fallback heuristic (or keep it as a
   floor for the case where OpenMP is unavailable).

PR 1 lands without any visible behaviour change; PR 2 lands the perf
win and is easy to revert if it misbehaves.

## Risks

- **OpenMP-vs-OpenBLAS oversubscription.** OpenBLAS already
  parallelises dgemm internally; running N threads each doing dgemm
  with the same OMP pool can thrash. Mitigation: set
  `OPENBLAS_NUM_THREADS=1` inside the parallel region (the existing
  `_fast_s` threaded path already does this with a context
  manager); inside Cython we would need an equivalent guard.

- **Per-thread bitgen seed quality.** If we naively use
  `seed + thread_id` instead of `SeedSequence.spawn`, adjacent
  streams can correlate. Always go through `SeedSequence`.

- **macOS LLVM OpenMP.** `libomp` on macOS sometimes needs an
  explicit dependency hint. `meson_options.txt` already marks OMP
  as best-effort. Verify the CI macOS jobs still build after the
  parallel patch lands.

## Estimate

PR 1 (refactor): 1 session.
PR 2 (parallel): 1-2 sessions including bench verification on the
synthetic n=10000, p=50 case and a fresh `docs/bench-report.md`.

The expected win is roughly N-way speedup on the resample loop at
large n, modulo BLAS-thread interaction. At n=5000, p=30 (current
engine_c is 424 ms), 4 threads should get us to ~120-150 ms which is
~0.12x R: a clear win over the current 0.42x.
