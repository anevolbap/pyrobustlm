# pylmrob

[![CI](https://github.com/anevolbap/pyrobustlm/actions/workflows/ci.yml/badge.svg)](https://github.com/anevolbap/pyrobustlm/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/pylmrob.svg)](https://pypi.org/project/pylmrob/)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](LICENSE)

Python port of the `lmrob` MM-estimator from R's
[`robustbase`](https://cran.r-project.org/package=robustbase) package.

Coefficient agreement with R on the classical validation datasets is within
`rtol=1e-3` for both estimates and covariance; see
[`docs/numerical-notes.md`](docs/numerical-notes.md) for known divergences
and [`docs/bench-report.md`](docs/bench-report.md) for wall-clock numbers.

Status: alpha. Public API is provisional.

## Install

```bash
pip install pylmrob
```

From source (requires Python 3.10+ and a C compiler):

```bash
git clone https://github.com/anevolbap/pyrobustlm
cd pyrobustlm
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
```

More: [`docs/quickstart.md`](docs/quickstart.md),
[`docs/examples/`](docs/examples/),
[`docs/porting-from-r.md`](docs/porting-from-r.md).

## Comparison with `statsmodels`

[`statsmodels.robust.RLM`](https://www.statsmodels.org/stable/rlm.html)
implements M-estimators (Huber, Tukey, Hampel) but not MM-estimators. If you
need MM (high breakdown + high efficiency), use `pylmrob`. If plain M is
enough, `RLM` has a more mature API.

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for the dev setup, test layout, and
PR guidelines.

## References

- Yohai, V. J. (1987). High Breakdown-Point and High Efficiency Robust
  Estimates for Regression. *Annals of Statistics*, 15(2).
  [doi:10.1214/aos/1176350366](https://doi.org/10.1214/aos/1176350366).
- Salibian-Barrera, M. and Yohai, V. J. (2006). A Fast Algorithm for
  S-Regression Estimates. *JCGS*, 15(2).
  [doi:10.1198/106186006X113629](https://doi.org/10.1198/106186006X113629).
- Koller, M. and Stahel, W. A. (2011). Sharpening Wald-type Inference in
  Robust Regression for Small Samples. *CSDA*, 55(8).
  [doi:10.1016/j.csda.2011.02.014](https://doi.org/10.1016/j.csda.2011.02.014).
- Koller, M. and Stahel, W. A. (2017). Nonsingular Subsampling for Regression
  S-estimators with Categorical Predictors. *Computational Statistics*,
  32(2). [doi:10.1007/s00180-016-0679-x](https://doi.org/10.1007/s00180-016-0679-x).
- Maronna, R. A. and Yohai, V. J. (2000). Robust regression with both
  continuous and categorical predictors. *JSPI*, 89(1-2).
  [doi:10.1016/S0378-3758(99)00208-6](https://doi.org/10.1016/S0378-3758(99)00208-6).

## License and citation

GPL-3.0-or-later, matching `robustbase`. See [`LICENSE`](LICENSE) and
[`NOTICE`](NOTICE). Citation metadata in [`CITATION.cff`](CITATION.cff)
(GitHub renders a "Cite this repository" button).
