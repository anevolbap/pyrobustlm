# SPDX-License-Identifier: GPL-3.0-or-later
"""Sphinx configuration for pylmrob."""

from __future__ import annotations

import os
import sys
from datetime import datetime, timezone
from pathlib import Path

# Headless matplotlib. The notebooks under docs/notebooks/ are executed
# at build time via myst-nb; we must pick the Agg backend BEFORE any
# notebook imports pyplot, otherwise CI / RTD pick the X11 backend and
# crash on no $DISPLAY.
try:
    import matplotlib as _mpl

    _mpl.use("Agg")
except ImportError:
    # matplotlib is in the ``docs`` optional-deps group; if a contributor
    # builds the docs from a minimal install we'd rather give a clear
    # sphinx error than fail here.
    pass

# Make the source tree importable so autodoc can introspect the package.
_repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo_root / "src"))

from pylmrob._version import __version__ as _pkg_version  # noqa: E402

# ---------------------------------------------------------------------------
# Project metadata
# ---------------------------------------------------------------------------
project = "pylmrob"
author = "anevolbap"
copyright = f"{datetime.now(timezone.utc).year}, {author}"

version = _pkg_version
release = _pkg_version

# ---------------------------------------------------------------------------
# General configuration
# ---------------------------------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    # myst_nb replaces myst_parser at runtime — it's a superset that
    # also handles ``{code-cell}`` blocks in markdown files. Listing
    # myst_parser too would cause a double-registration error.
    "myst_nb",
]

# myst-nb: execute notebooks at build time. ``cache`` means unchanged
# notebooks reuse last build's output; only first build (or after the
# notebook changes) actually runs the cells.
nb_execution_mode = "cache"
# Slow Azure runners can push notebook 03 (contamination sweep, 11 x 25
# lmrob fits) past the 120s default; 300s leaves room without masking a
# genuinely runaway cell.
nb_execution_timeout = 300
nb_execution_raise_on_error = True
# Notebooks under docs/notebooks/ should always be MyST markdown
# (no .ipynb); the next line is a belt-and-braces declaration.
nb_render_image_options = {"width": "100%"}

# myst_nb registers .md and .ipynb itself; .rst is restructuredtext.
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
}

master_doc = "index"

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "api/_autosummary/_*",
    # Legacy long-form notes that live as plain markdown on disk and on
    # GitHub; not part of the rendered site (yet).
    "research-notes.md",
    "r-source-map.md",
    "bench-report.md",
]

# Autodoc
autodoc_default_options = {
    "members": True,
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"
autosummary_generate = True

# Napoleon (numpy-style docstrings).
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_param = True
napoleon_use_rtype = True

# Intersphinx
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
}

# MyST: enable a few useful extensions.
#
# ``dollarmath`` is what hands ``$...$`` and ``$$...$$`` blocks to
# mathjax; without it the math in theory.md renders as literal text.
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "smartquotes",
]

nitpicky = False  # keep warnings actionable; flip to True before publish.

# ---------------------------------------------------------------------------
# HTML output
# ---------------------------------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_show_sourcelink = False

html_theme_options = {
    "navigation_depth": 3,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "style_external_links": True,
}

# Suppress the noisy warning about missing _static dir on a fresh checkout.
if not (Path(__file__).parent / "_static").exists():
    os.makedirs(Path(__file__).parent / "_static", exist_ok=True)
