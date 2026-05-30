# SPDX-License-Identifier: GPL-3.0-or-later
"""pylmrob: Python port of robustbase::lmrob (MM regression).

The supported public surface is listed in ``__all__`` below. It follows
the deprecation policy in ``docs/policy.md``: changes go through one
minor release of ``DeprecationWarning`` before removal.

Anything imported from a leading-underscore submodule (``pylmrob._*``)
is internal and may change in any release.
"""

from __future__ import annotations

from pylmrob._fast_s import make_generator as make_generator
from pylmrob._version import __version__
from pylmrob.anova import anova
from pylmrob.bootstrap import bootstrap
from pylmrob.control import Control
from pylmrob.lmrob import LmRob, lmrob
from pylmrob.results import LmRobResults
from pylmrob.rng import RState as RState
from pylmrob.rng import r_norm_rand as r_norm_rand
from pylmrob.rng import r_qnorm as r_qnorm
from pylmrob.rng import r_sample_noreplace as r_sample_noreplace
from pylmrob.rng import r_set_seed as r_set_seed
from pylmrob.rng import r_subsample_nonsingular as r_subsample_nonsingular

# Stable public surface. The ``r_*`` helpers are kept importable from
# the top level for backwards compatibility but are not in ``__all__``
# because they live in the advanced ``pylmrob.rng`` namespace and are
# meant for the R-bridge validation workflow rather than day-to-day use.
__all__ = [
    "Control",
    "LmRob",
    "LmRobResults",
    "__version__",
    "anova",
    "bootstrap",
    "lmrob",
]
