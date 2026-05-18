# SPDX-License-Identifier: GPL-3.0-or-later
"""Phase 0 smoke tests.

These confirm the package, the version string, and the compiled stub
extension all import. They will keep guarding regressions of the build
system as later phases land.
"""

from __future__ import annotations

import importlib

import numpy as np

import pylmrob

_stub = importlib.import_module("pylmrob._core._stub")


def test_version_is_string() -> None:
    assert isinstance(pylmrob.__version__, str)
    assert pylmrob.__version__


def test_public_names_present() -> None:
    assert hasattr(pylmrob, "lmrob")
    assert hasattr(pylmrob, "Control")
    assert hasattr(pylmrob, "LmRob")


def test_stub_extension_loads() -> None:
    assert _stub.hello() == "pylmrob._core._stub OK"


def test_stub_vec_norm() -> None:
    x = np.arange(1.0, 5.0)
    expected = float(np.linalg.norm(x))
    assert abs(_stub.vec_norm(x) - expected) < 1e-12


def test_control_preset_round_trip() -> None:
    ctrl = pylmrob.Control.preset("KS2014", nResample=100)
    assert ctrl.setting == "KS2014"
    assert ctrl.nResample == 100
