# SPDX-License-Identifier: GPL-3.0-or-later
"""pyrobustlm: Python port of robustbase::lmrob (MM regression).

See plan.md for the full roadmap. The public surface listed below is
provisional; most symbols are stubs until their phase lands.
"""

from __future__ import annotations

from pyrobustlm._version import __version__
from pyrobustlm.control import Control
from pyrobustlm.lmrob import LmRob, lmrob

__all__ = [
    "Control",
    "LmRob",
    "__version__",
    "lmrob",
]
