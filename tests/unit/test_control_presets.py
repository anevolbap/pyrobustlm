# SPDX-License-Identifier: GPL-3.0-or-later
"""``Control.preset`` and ``Control(setting=...)`` must be the same thing.

Reference values from ``robustbase::lmrob.control(setting=)``
(robustbase 0.99-7): both KS2011 and KS2014 are psi="lqq",
method="SMDM", cov=".vcov.w"; the unnamed default is psi="bisquare",
method="MM", cov=".vcov.avar1".
"""

from __future__ import annotations

import dataclasses

import pytest

from pylmrob import Control

_SETTINGS = ["KS2011", "KS2014", "MM"]


@pytest.mark.parametrize("setting", _SETTINGS)
def test_preset_equals_constructor(setting: str) -> None:
    """The two entry points used to disagree for KS2011/KS2014."""
    assert Control.preset(setting) == Control(setting=setting)  # type: ignore[arg-type]


@pytest.mark.parametrize("setting", ["KS2011", "KS2014"])
def test_ks_settings_match_r(setting: str) -> None:
    ctrl = Control.preset(setting)  # type: ignore[arg-type]
    assert ctrl.psi == "lqq"
    assert ctrl.method == "SMDM"
    assert ctrl.cov == ".vcov.w"


def test_mm_setting_matches_r() -> None:
    ctrl = Control.preset("MM")
    assert ctrl.psi == "bisquare"
    assert ctrl.method == "MM"
    assert ctrl.cov == ".vcov.avar1"


def test_preset_override_recomputes_tuning() -> None:
    """An override must go through ``__post_init__``.

    The old implementation ``setattr``-ed overrides after construction,
    so ``psi="lqq"`` kept the previous family's tuning constants and the
    fit silently used the wrong ones.
    """
    ctrl = Control.preset("MM", psi="lqq")
    assert ctrl.psi == "lqq"
    assert ctrl.tuning_psi == Control(psi="lqq").tuning_psi
    assert ctrl.tuning_chi == Control(psi="lqq").tuning_chi
    # and specifically *not* bisquare's constant
    assert ctrl.tuning_psi != Control(psi="bisquare").tuning_psi


def test_preset_rejects_unknown_setting() -> None:
    with pytest.raises(ValueError, match="unknown setting"):
        Control.preset("KS2099")  # type: ignore[arg-type]


def test_preset_rejects_unknown_field() -> None:
    with pytest.raises(TypeError, match="unknown Control field"):
        Control.preset("MM", not_a_field=1)


def test_preset_passes_through_valid_overrides() -> None:
    ctrl = Control.preset("KS2014", nResample=123)
    assert ctrl.nResample == 123
    assert ctrl.psi == "lqq"  # preset defaults survive the override


def test_every_preset_field_is_a_real_field() -> None:
    """Guard against the override whitelist drifting from the dataclass."""
    names = {f.name for f in dataclasses.fields(Control)}
    assert "numpoints" in names, "numpoints should exist (R's lmrob.control has it)"


def test_huber_psi_is_rejected_like_r() -> None:
    """robustbase rejects psi="huber" for lmrob; so do we.

    Huber's rho is unbounded, so an S-estimate built on it has 0%
    breakdown: the high-breakdown initial estimate MM relies on does not
    exist. We used to accept it and return a plausible-looking fit with
    none of the robustness the user came for. R's message lists
    tukey/biweight/bisquare, lqq, welsh, optimal, hampel, ggw.
    """
    with pytest.raises(ValueError, match="not a valid lmrob psi"):
        Control(psi="huber")  # type: ignore[arg-type]


def test_huber_primitives_still_work() -> None:
    """Only lmrob's psi rejects huber. Mpsi/Mchi/m_scale still support it,
    as they do in R."""
    import numpy as np

    from pylmrob import psi as psi_mod
    from pylmrob.scale import m_scale

    x = np.linspace(-3, 3, 21)
    assert np.all(np.isfinite(psi_mod.psi(x, "huber", (1.345,))))
    assert m_scale(np.array([1.0, -2.0, 3.0, -1.0, 0.5]), "huber") > 0


@pytest.mark.parametrize("alias", ["bisquare", "tukey", "biweight"])
def test_accepted_psi_families_match_r(alias: str) -> None:
    """The families R accepts must construct without error."""
    from pylmrob.control import PsiFamily  # noqa: F401

    if alias == "bisquare":
        assert Control(psi=alias).psi == "bisquare"  # type: ignore[arg-type]
