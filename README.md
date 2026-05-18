# pylmrob

[![CI](https://github.com/anevolbap/pylmrob/actions/workflows/ci.yml/badge.svg)](https://github.com/anevolbap/pylmrob/actions/workflows/ci.yml)
[![Docs](https://readthedocs.org/projects/pylmrob/badge/?version=latest)](https://pylmrob.readthedocs.io/)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](LICENSE)

Python port of the `lmrob` MM-estimator from R's
[`robustbase`](https://cran.r-project.org/package=robustbase) package.

Documentation: [pylmrob.readthedocs.io](https://pylmrob.readthedocs.io/).

Status: v0.5.11 (alpha). End-to-end pipeline matches R's `lmrob` element-wise
on the classical robust regression datasets (`stackloss`, `coleman`,
`delivery`, `aircraft`, `phosphor`, etc.) within `rtol=1e-3` for both
coefficients and covariance. Median wall-clock is **0.93x R** across the
34-case bench corpus. See [`plan.md`](plan.md) for the full roadmap and
[`docs/numerical-notes.md`](docs/numerical-notes.md) for known
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
- Per-case `weights` argument: matches R element-wise on the NumPy
  path. Non-trivial weights force ``Control(engine_c=False)`` for now.
- Default `Control()` runs the monolithic Cython engine (single nogil
  C block for fast-S, MM, D-scale, and vcov); see
  [Performance](#performance).

## What does not work yet

- Bit-identical reproducibility with R's MT RNG. We use NumPy's PCG64
  (waived in `plan.md` §5.2; coefficients agree with R within
  basin-of-attraction tolerances documented in
  [`docs/numerical-notes.md`](docs/numerical-notes.md)).

## Performance

Median wall-clock vs R across the 34-case bench corpus is **0.93x**;
pylmrob is faster than R on more than half of cases (down to 0.32x R
on `n>=2000` fits) and within 1-2x R on small classical datasets.

The default `Control()` runs the fit through a monolithic Cython engine
(fast-S, MM, D-scale, and vcov in one nogil C block). At larger n
(`n*p^2 >= 100k`) `lmrob()` auto-falls-back to a threaded NumPy path
because the monolithic kernel is a single non-parallel call. To force
the legacy NumPy path explicitly, pass `Control(engine_c=False)`.

See [`docs/bench-report.md`](docs/bench-report.md) for the full per-case
table and [`docs/engine_c.md`](docs/engine_c.md) for the trade-off
discussion.

## Install (from source)

Requires Python 3.10+, a C compiler, and (for the validation harness) a
local R install with `robustbase`.

```bash
git clone https://github.com/anevolbap/pylmrob
cd pylmrob
uv venv
uv pip install --no-build-isolation -e ".[dev]"
pytest
```

## Quickstart

```python
import pandas as pd
from pylmrob import lmrob, Control

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
[pylmrob.readthedocs.io](https://pylmrob.readthedocs.io/) and
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
`lmrob` provides. If you need MM, use `pylmrob`. If you only need M
with no contamination concerns, `RLM` is fine and has a more mature API.

## License

GPL-3.0-or-later, matching `robustbase`. See [`LICENSE`](LICENSE) and
[`NOTICE`](NOTICE).

## Citing

If you use `pylmrob` in research, please cite the project alongside
`robustbase`. Citation metadata is in [`CITATION.cff`](CITATION.cff)
(GitHub shows a "Cite this repository" button on the sidebar). The
references include the key papers: Yohai (1987) on MM-estimators,
Salibian-Barrera & Yohai (2006) on the fast-S algorithm, and
Koller & Stahel (2011) on the small-sample Wald correction implemented
in `vcov.w`.
