# SPDX-License-Identifier: GPL-3.0-or-later
"""Every tuning constant we hardcode, read back from R.

Three bugs in a row were constants transcribed by hand from R output at
reduced precision:

* the D-step ``kappa`` table, whose ggw case-4 entry was a copy of case
  1, putting the D-scale 3.8% off on every ggw fit;
* ``_DEFAULT_TUNING_CHI["bisquare"]``, which was ``1.547645`` against
  R's ``1.54764`` -- a 3.231e-06 relative error that showed up as the
  "scale relative error: median 3.23e-06" floor in every revision of
  ``docs/bench-report.md``, and was misdiagnosed as irreducible LAPACK
  noise for several releases;
* the lqq psi mid-constant, ``0.9826779`` against ``0.9822707``, which
  silently disabled a precomputed correction factor in ``inference.py``
  because its ``np.allclose`` guard stopped matching.

None was caught by a test, because every test that used a constant
imported the same wrong constant it was checking. This file breaks that
loop: it asks R and compares.

The source of truth is ``robustbase:::.Mpsi.tuning.default`` /
``.Mchi.tuning.default`` (what ``control.py`` cites) pushed through
``.psi.conv.cc``, which is identity for the plain families, maps lqq to
its ``(b, c, s)`` triple, and maps ggw to its SET_ABC_GGW case index.

Gated on a working ``Rscript`` with ``robustbase``; CI installs both on
the Linux job.
"""

from __future__ import annotations

import json
import subprocess

import numpy as np
import pytest

from pylmrob.control import _DEFAULT_TUNING_CHI, _DEFAULT_TUNING_PSI

# R has no chi default for huber: its rho is unbounded, so there is no
# 50%-breakdown constant to define. Our ``_DEFAULT_TUNING_CHI["huber"]``
# is therefore ours alone, and huber is not a supported S/MM chi family
# (it is absent from the Cython family tables). psi is still checked.
_NO_R_CHI = {"huber"}

_FAMILIES = sorted(set(_DEFAULT_TUNING_PSI) & set(_DEFAULT_TUNING_CHI))
_CHI_FAMILIES = [f for f in _FAMILIES if f not in _NO_R_CHI]


def _r_has_robustbase() -> bool:
    try:
        proc = subprocess.run(
            ["Rscript", "-e", "library(robustbase)"],
            capture_output=True,
            text=True,
            timeout=60,
        )
        return proc.returncode == 0
    except Exception:
        return False


pytestmark = pytest.mark.skipif(
    not _r_has_robustbase(),
    reason="Rscript + robustbase not available",
)

_R_SCRIPT = """
suppressMessages({library(robustbase); library(jsonlite)})
fams <- commandArgs(trailingOnly = TRUE)
out <- list()
for (ps in fams) {
  one <- list()
  tp <- tryCatch(robustbase:::.Mpsi.tuning.default(ps), error = function(e) NULL)
  tc <- tryCatch(robustbase:::.Mchi.tuning.default(ps), error = function(e) NULL)
  if (!is.null(tp)) one$psi <- as.numeric(robustbase:::.psi.conv.cc(ps, tp))
  if (!is.null(tc)) one$chi <- as.numeric(robustbase:::.psi.conv.cc(ps, tc))
  out[[ps]] <- one
}
cat(toJSON(out, digits = 17, auto_unbox = FALSE))
"""


def _r_tuning() -> dict[str, dict[str, list[float]]]:
    proc = subprocess.run(
        ["Rscript", "-e", _R_SCRIPT, *_FAMILIES],
        capture_output=True,
        text=True,
        timeout=180,
    )
    if proc.returncode != 0:  # pragma: no cover - surfaced as a failure
        pytest.fail(f"Rscript failed:\n{proc.stderr}")
    return json.loads(proc.stdout)


@pytest.fixture(scope="module")
def r_tuning() -> dict[str, dict[str, list[float]]]:
    return _r_tuning()


def _compare(family: str, kind: str, ours_raw: object, r_tuning: dict) -> None:
    entry = r_tuning.get(family, {})
    if kind not in entry:
        pytest.fail(f"{family}: R defines no {kind} default; update _NO_R_CHI if intended")
    ours = np.atleast_1d(np.asarray(ours_raw, dtype=float)).ravel()
    theirs = np.atleast_1d(np.asarray(entry[kind], dtype=float)).ravel()
    assert ours.shape == theirs.shape, f"{family} {kind}: length {ours.shape} vs R {theirs.shape}"
    np.testing.assert_allclose(
        ours, theirs, rtol=1e-15, atol=0.0, err_msg=f"{family}: tuning_{kind} disagrees with R"
    )


@pytest.mark.parametrize("family", _FAMILIES)
def test_tuning_psi_matches_r(family: str, r_tuning: dict) -> None:
    _compare(family, "psi", _DEFAULT_TUNING_PSI[family], r_tuning)


@pytest.mark.parametrize("family", _CHI_FAMILIES)
def test_tuning_chi_matches_r(family: str, r_tuning: dict) -> None:
    """The assertion that would have caught the bisquare ``1.547645`` typo."""
    _compare(family, "chi", _DEFAULT_TUNING_CHI[family], r_tuning)


def test_r_really_has_no_huber_chi(r_tuning: dict) -> None:
    """Pin the reason huber is exempt, so the exemption cannot rot."""
    for family in _NO_R_CHI:
        assert "chi" not in r_tuning.get(family, {}), (
            f"R now defines a chi default for {family}; drop it from _NO_R_CHI "
            "and check our value against it"
        )


def test_duplicate_tables_agree() -> None:
    """The same constants live in several modules; they must not drift.

    ``inference.py`` carried two lqq tuning tables that disagreed with
    each other and with ``Control``, and the stale one silently disabled
    a precomputed correction factor.
    """
    from pylmrob import psi as psi_mod
    from pylmrob.scale import _DEFAULT_K_CHI

    for family in _FAMILIES:
        expected = np.atleast_1d(np.asarray(_DEFAULT_TUNING_CHI[family], dtype=float)).ravel()
        for label, table in (
            ("scale._DEFAULT_K_CHI", _DEFAULT_K_CHI),
            ("psi._PSI_TUNING_DEFAULT_CHI", getattr(psi_mod, "_PSI_TUNING_DEFAULT_CHI", {})),
        ):
            if family not in table:
                continue
            got = np.atleast_1d(np.asarray(table[family], dtype=float)).ravel()
            if got.shape != expected.shape:
                pytest.fail(
                    f"{family}: {label} has length {got.shape}, Control has {expected.shape}"
                )
            np.testing.assert_allclose(
                got,
                expected,
                rtol=1e-15,
                atol=0.0,
                err_msg=f"{family}: {label} drifted from Control._DEFAULT_TUNING_CHI",
            )
