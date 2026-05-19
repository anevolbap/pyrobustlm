# Contributing

Thanks for your interest. `pylmrob` is alpha; the public API is provisional
and may change without notice. Bug reports, R-to-Python divergence reports,
and small fixes are all welcome.

## Dev setup

Requires Python 3.10+, a C compiler, `uv`, and (optionally) a local R with
`robustbase` for the validation harness.

```bash
git clone https://github.com/anevolbap/pyrobustlm
cd pyrobustlm
uv venv
uv pip install --no-build-isolation -e ".[dev]"
pre-commit install
```

For the R-comparison tests:

```bash
uv pip install --no-build-isolation -e ".[dev,validation]"
```

For the docs build:

```bash
uv pip install --no-build-isolation -e ".[docs]"
make -C docs html
```

## Tests

The suite is in `tests/`, split into `unit/`, `integration/`, `property/`
(Hypothesis), `reference/` (per-function R comparisons), and `validation/`
(end-to-end R comparisons; requires `rpy2`).

```bash
pytest                              # full suite
pytest tests/unit -x                # one layer, fail fast
pytest -k "stackloss"               # filter by name
pytest --no-cov                     # skip coverage when iterating
```

## Style

`ruff` and `pre-commit` enforce formatting; the hooks run on commit.
Conventional-commits prefixes for messages (`feat:`, `fix:`, `refactor:`,
`chore:`, `docs:`).

## Issues and PRs

Issues: <https://github.com/anevolbap/pyrobustlm/issues>. Please include
`pylmrob.__version__`, Python version, OS, and a minimal reproducer. If the
bug is a numerical disagreement with R, include the R call and output.

PRs: keep them small and focused. Tests should pass locally before pushing;
CI runs the full matrix on push.

## Reference reading

Internal guides:

- [`docs/theory.md`](docs/theory.md) - M / S / MM estimators.
- [`docs/porting-from-r.md`](docs/porting-from-r.md) - mapping from
  `robustbase::lmrob` API to `pylmrob`.
- [`docs/r-source-map.md`](docs/r-source-map.md) - which Python file
  ports which R / C file in `robustbase`.
- [`docs/numerical-notes.md`](docs/numerical-notes.md) - known divergences
  from R and their causes.

Papers (DOI links):

- Yohai, V. J. (1987). High Breakdown-Point and High Efficiency Robust
  Estimates for Regression. *Annals of Statistics*, 15(2).
  [doi:10.1214/aos/1176350366](https://doi.org/10.1214/aos/1176350366).
- Salibian-Barrera, M. and Yohai, V. J. (2006). A Fast Algorithm for
  S-Regression Estimates. *JCGS*, 15(2).
  [doi:10.1198/106186006X113629](https://doi.org/10.1198/106186006X113629).
- Koller, M. and Stahel, W. A. (2011). Sharpening Wald-type Inference in
  Robust Regression for Small Samples. *CSDA*, 55(8).
  [doi:10.1016/j.csda.2011.02.014](https://doi.org/10.1016/j.csda.2011.02.014).
- Koller, M. and Stahel, W. A. (2017). Nonsingular Subsampling for
  Regression S-estimators with Categorical Predictors. *Computational
  Statistics*, 32(2).
  [doi:10.1007/s00180-016-0679-x](https://doi.org/10.1007/s00180-016-0679-x).
- Maronna, R. A. and Yohai, V. J. (2000). Robust regression with both
  continuous and categorical predictors. *JSPI*, 89(1-2).
  [doi:10.1016/S0378-3758(99)00208-6](https://doi.org/10.1016/S0378-3758(99)00208-6).
