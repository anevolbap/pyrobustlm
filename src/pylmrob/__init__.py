# SPDX-License-Identifier: GPL-3.0-or-later
"""pylmrob: Python port of robustbase::lmrob (MM regression).

See plan.md for the full roadmap. The public surface listed below is
provisional; most symbols are stubs until their phase lands.
"""

from __future__ import annotations

from pylmrob._fast_s import make_generator
from pylmrob._version import __version__
from pylmrob.anova import anova
from pylmrob.bootstrap import bootstrap
from pylmrob.control import Control
from pylmrob.lmrob import LmRob, lmrob

__all__ = [
    "Control",
    "LmRob",
    "__version__",
    "anova",
    "bootstrap",
    "lmrob",
    "make_generator",
]
