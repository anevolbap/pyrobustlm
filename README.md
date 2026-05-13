# pyrobustlm

[![CI](https://github.com/anevolbap/pyrobustlm/actions/workflows/ci.yml/badge.svg)](https://github.com/anevolbap/pyrobustlm/actions/workflows/ci.yml)
[![Docs](https://readthedocs.org/projects/pyrobustlm/badge/?version=latest)](https://pyrobustlm.readthedocs.io/)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](LICENSE)

Python port of the `lmrob` MM-estimator from R's
[`robustbase`](https://cran.r-project.org/package=robustbase) package.

Documentation: [pyrobustlm.readthedocs.io](https://pyrobustlm.readthedocs.io/).

Status: v0.5.5 (alpha). End-to-end pipeline matches R's `lmrob` element-wise
on the classical robust regression datasets (`stackloss`, `coleman`,
`delivery`, `aircraft`, `phosphor`, etc.) within `rtol=1e-3` for both
coefficients and covariance. See [`plan.md`](plan.md) for the full roadmap
and [`docs/numerical-notes.md`](docs/numerical-notes.md) for known
divergences from R.

## What works

- All six psi families (`bisquare`, `huber`, `hampel`, `optimal`, `lqq`,
  `ggw`) match R's `Mpsi`/`Mchi`/`Mwgt` to `rtol=1e-12` (1e-7 for ggw rho's
  polynomial approximation).
- M-scale matches R's `find_scale` to `rtol=1e-9`.
- Fast-S resampling + MM iteration produces final coefficients matching R
  on the validation corpus (10 classical datasets) within per-dataset
  tolerances ranging `1e-3` (well-conditioned) to `1e-1` (small-n).
- Sandwich covariance estimator (`vcov_avar1`).
- `Control` presets matching R's `lmrob.control(setting={"KS2014","KS2011","MM"})`.
- Affine equivariance, scale equivariance, regression equivariance verified
  via Hypothesis property tests.
- `init="M-S"` for factor designs (Maronna-Yohai 2000, 4-phase port of
  robustbase's `R_lmrob_M_S`).
- D-scale refinement (Koller & Stahel 2014) for `setting="KS2014"` and
  `setting="KS2011"`; full SMDM pipeline matches R element-wise.
- `vcov.w` with all five `cov.corrfact` branches plus the Huber
  finite-sample correction, matching R within `rtol=1e-3`.
- `summary()` (coefficient table, robust R-squared) and `anova()`
  (Wald + Deviance) matching R element-wise.
- `Control(engine_c=True)` opt-in monolithic Cython engine; see
  [Performance](#performance).

## What does not work yet

- Per-case `weights` argument (raises `NotImplementedError`).
- Bit-identical reproducibility with R's MT RNG. We use NumPy's PCG64
  (waived in `plan.md` §5.2; coefficients agree with R within
  basin-of-attraction tolerances documented in
  [`docs/numerical-notes.md`](docs/numerical-notes.md)).

## Performance

The default path is pure-NumPy + Cython kernels for the inner loops and
runs about 5-15x R wall-clock on small-n problems (see
[`docs/bench-report.md`](docs/bench-report.md)).

`Control(engine_c=True)` opts into a monolithic Cython lmrob engine that
runs the entire fit (fast-S, MM, D-scale, vcov) in one nogil C block.
On a single-threaded BLAS box, phosphor MM fits in ~3.6 ms vs R's
~3 ms (about 1.2x R); KS settings 6-10 ms vs R's 4-7 ms. See
[`docs/engine_c.md`](docs/engine_c.md) for the trade-off and the
reproducible bench script `scripts/bench_engine_c.py`.

## Install (from source)

Requires Python 3.10+, a C compiler, and (for the validation harness) a
local R install with `robustbase`.

```bash
git clone https://github.com/anevolbap/pyrobustlm
cd pyrobustlm
uv venv
uv pip install --no-build-isolation -e ".[dev]"
pytest
```

## Quickstart

```python
import pandas as pd
from pyrobustlm import lmrob, Control

df = pd.read_csv("stackloss.csv")
fit = lmrob(
    "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
    df,
    control=Control(setting="KS2014"),
    seed=42,
)
print(fit.summary())

# Coefficients match R within rtol=1e-3:
# (Intercept) -41.5246
# Air.Flow      0.9388
# Water.Temp    0.5796
# Acid.Conc.   -0.1129
# scale = 1.9123

# Predict on new rows; the formula spec is re-applied automatically.
new_rows = pd.DataFrame({
    "Air.Flow":    [60, 80],
    "Water.Temp":  [20, 25],
    "Acid.Conc.":  [85, 90],
})
fit.predict(new_rows)

# Confidence intervals:
fit.confint(level=0.95)
```

## Documentation

The Sphinx docs are hosted at
[pyrobustlm.readthedocs.io](https://pyrobustlm.readthedocs.io/) and
built on every push (config in [`.readthedocs.yaml`](.readthedocs.yaml)).
To build locally:

```bash
uv pip install -e ".[docs]" --no-build-isolation
cd docs && uv run sphinx-build -b html -W . _build/html
open _build/html/index.html
```

Pages: quickstart, API reference, R-to-Python porting cheatsheet, and the
numerical-notes log of documented divergences from R.

## Comparison with `statsmodels`

[`statsmodels.robust.RLM`](https://www.statsmodels.org/stable/rlm.html)
implements M-estimators (Huber, Tukey, Hampel) but **not** MM-estimators.
It will not give you the high-breakdown-point + high-efficiency combo that
`lmrob` provides. If you need MM, use `pyrobustlm`. If you only need M
with no contamination concerns, `RLM` is fine and has a more mature API.

## License

GPL-3.0-or-later, matching `robustbase`. See [`LICENSE`](LICENSE) and
[`NOTICE`](NOTICE).

## Citing

If you use `pyrobustlm` in research, please cite both `robustbase`'s
canonical reference (Maechler et al., 2024) and this project's GitHub URL
until a Zenodo DOI is minted.
