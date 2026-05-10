# SPDX-License-Identifier: GPL-3.0-or-later
"""Sphinx configuration for pyrobustlm."""

from __future__ import annotations

import os
import sys
from datetime import datetime, timezone
from pathlib import Path

# Make the source tree importable so autodoc can introspect the package.
_repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_repo_root / "src"))

from pyrobustlm._version import __version__ as _pkg_version  # noqa: E402

# ---------------------------------------------------------------------------
# Project metadata
# ---------------------------------------------------------------------------
project = "pyrobustlm"
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
    "myst_parser",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
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
myst_enable_extensions = [
    "colon_fence",
    "deflist",
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
