# SPDX-License-Identifier: GPL-3.0-or-later
"""Exact fits must return, not raise.

Noiseless data (``y = 1 + 2x``) drives the S-scale to zero, which is a
legitimate S solution: the fit passes through every point. Every
residual is then zero, so ``vcov_avar1``'s denominator
``mean(chi'(r0/s) * r0/s)`` is zero because ``chi'(0) * 0 == 0``, not
because anything went wrong.

We used to raise ``FloatingPointError`` there. R warns
("S-estimated scale == 0: Probably exact fit; check your data") and
returns coefficients with a zero covariance, so a user fitting a
noiseless design gets an answer. Verified against robustbase 0.99-7:

    coef: 1.000000 2.000000   scale: 0   vcov diag: 0, 0
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from pylmrob import Control, lmrob


def _noiseless() -> pd.DataFrame:
    x = np.arange(20, dtype=float)
    return pd.DataFrame({"x": x, "y": 1.0 + 2.0 * x})


def test_exact_fit_returns_instead_of_raising() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = lmrob("y ~ x", _noiseless(), control=Control(nResample=100), seed=0)

    np.testing.assert_allclose(fit.coef_, [1.0, 2.0], rtol=1e-9, atol=1e-9)
    assert fit.scale_ == 0.0


def test_exact_fit_covariance_is_zero_like_r() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = lmrob("y ~ x", _noiseless(), control=Control(nResample=100), seed=0)

    assert fit.cov_.shape == (2, 2)
    np.testing.assert_allclose(fit.cov_, np.zeros((2, 2)), atol=0.0)


def test_exact_fit_warns() -> None:
    """Silence would be worse than the old raise: a zero covariance means
    every standard error is zero, which a caller must not take at face
    value."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        lmrob("y ~ x", _noiseless(), control=Control(nResample=100), seed=0)

    messages = " ".join(str(w.message) for w in caught)
    assert "exact fit" in messages.lower(), messages


def test_ordinary_fit_is_unaffected() -> None:
    """The guard must not swallow a genuinely degenerate covariance."""
    rng = np.random.default_rng(0)
    n = 40
    x = rng.standard_normal(n)
    df = pd.DataFrame({"x": x, "y": 1.0 + 2.0 * x + rng.standard_normal(n) * 0.5})
    fit = lmrob("y ~ x", df, control=Control(nResample=200), seed=1)

    assert fit.scale_ > 0.0
    assert np.all(np.diag(fit.cov_) > 0.0)
