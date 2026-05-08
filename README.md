# pyrobustlm

Python port of the `lmrob` MM-estimator from R's
[`robustbase`](https://cran.r-project.org/package=robustbase) package.

Status: pre-alpha. See [`plan.md`](plan.md) for the full roadmap.

## Goals

- Feature parity with `robustbase::lmrob`: all settings (`KS2011`, `KS2014`,
  default MM), all init strategies (`S`, `M-S`, `L1`), all six psi families
  (`bisquare`, `huber`, `hampel`, `optimal`, `ggw`, `lqq`), and all covariance
  estimators.
- Numerical agreement with R within documented tolerances on a fixed
  validation corpus.
- Performance at parity with R on small problems and faster on large ones
  when multiple cores are available.

## Install (once published)

```bash
pip install pyrobustlm
```

## Quickstart

```python
import pandas as pd
from pyrobustlm import lmrob

df = pd.read_csv("stackloss.csv")
fit = lmrob("stack_loss ~ Air_Flow + Water_Temp + Acid_Conc", df)
print(fit.summary())
```

## Development

Requires Python 3.10+, a C compiler, and (for the validation harness) a local
R install with `robustbase`.

```bash
uv venv
uv pip install -e ".[dev]"
pytest
```

## License

GPL-3.0-or-later. See [`LICENSE`](LICENSE) and [`NOTICE`](NOTICE).
