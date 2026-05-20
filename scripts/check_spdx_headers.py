#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Pre-commit hook: every source file must carry an SPDX license header.

Pylmrob is GPL-3.0-or-later. Each file gets:

    # SPDX-License-Identifier: GPL-3.0-or-later

within the first 5 lines. Cython ``# cython: ...`` lines and shebangs
are allowed before it.
"""

from __future__ import annotations

import sys
from pathlib import Path

NEEDLE = "SPDX-License-Identifier: GPL-3.0-or-later"


def has_header(path: Path) -> bool:
    try:
        with path.open("r", encoding="utf-8") as fh:
            head = "".join(fh.readline() for _ in range(5))
    except OSError:
        return False
    return NEEDLE in head


def main(argv: list[str]) -> int:
    missing = [p for p in argv[1:] if not has_header(Path(p))]
    if missing:
        sys.stderr.write(
            f"\nMissing SPDX header ({NEEDLE!r}) in {len(missing)} file(s):\n"
        )
        for p in missing:
            sys.stderr.write(f"  {p}\n")
        sys.stderr.write(
            "\nAdd the line within the first 5 lines of the file. For Python /"
            " Cython:\n\n    # SPDX-License-Identifier: GPL-3.0-or-later\n\n"
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
