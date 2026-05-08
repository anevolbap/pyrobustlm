# SPDX-License-Identifier: GPL-3.0-or-later
"""Diagnostic plots and influence statistics for ``lmrob`` fits.

Phase 9.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pyrobustlm.results import LmRobResults


def plot(results: LmRobResults) -> object:
    raise NotImplementedError("diagnostics.plot — Phase 9")


def cooks_distance(results: LmRobResults, robust: bool = True) -> object:
    raise NotImplementedError("diagnostics.cooks_distance — Phase 9")


def hatvalues(results: LmRobResults) -> object:
    raise NotImplementedError("diagnostics.hatvalues — Phase 9")
