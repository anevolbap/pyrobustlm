#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Pre-commit hook: ``meson.build`` ``version`` and ``_version.py``
``__version__`` must agree.

This is the bug that took out the v0.5.16 and v0.5.17 publishes: the
Python version was bumped, but ``meson.build`` still said the old
value, so cibuildwheel produced wheels named for the stale version and
PyPI rejected them with "File already exists" for the v0.5.15 file.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
VERSION_PY = REPO / "src" / "pylmrob" / "_version.py"
MESON = REPO / "meson.build"

_PY_RE = re.compile(r'__version__\s*=\s*"([^"]+)"')
_MESON_RE = re.compile(r"version:\s*'([^']+)'")


def _read(path: Path, regex: re.Pattern[str], label: str) -> str:
    m = regex.search(path.read_text())
    if m is None:
        sys.stderr.write(f"could not find version in {label}\n")
        sys.exit(2)
    return m.group(1)


def main() -> int:
    py = _read(VERSION_PY, _PY_RE, "_version.py")
    meson = _read(MESON, _MESON_RE, "meson.build")
    if py != meson:
        sys.stderr.write(
            f"\nVersion mismatch:\n"
            f"  src/pylmrob/_version.py __version__ = {py!r}\n"
            f"  meson.build version           = {meson!r}\n\n"
            f"Both must be the same. Bump them together.\n"
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
