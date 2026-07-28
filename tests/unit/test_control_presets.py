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
