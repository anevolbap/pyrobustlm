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

# Compiled extensions ``_stub`` and ``_psi`` are loaded lazily by their
# callers via ``importlib`` so that static type checkers don't trip on the
# Cython-built modules at analysis time.
