# SPDX-License-Identifier: GPL-3.0-or-later
"""Compiled kernels.

The submodules in this package are Cython extensions and are populated
phase by phase:

- ``_psi``           Phase 2
- ``_scale``         Phase 3
- ``_fast_s``        Phase 4
- ``_mm``            Phase 6
- ``_lapack_helpers``  Phase 4
- ``_ctables``       Phase 2

For Phase 0 only ``_stub`` exists; it exists to verify the build toolchain.
"""

from pyrobustlm._core import _stub  # noqa: F401
