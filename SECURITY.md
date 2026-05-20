# Security Policy

## Supported versions

`pylmrob` is in alpha (0.x). Only the latest released minor version on PyPI
receives security fixes; older versions are not patched.

## Reporting a vulnerability

If you find a security issue, please **do not** open a public issue. Instead:

1. Email the maintainer at <vena.pablo@gmail.com> with the subject `pylmrob security`.
2. Include a minimal reproducer (data + the exact call), the version affected,
   and your assessment of the impact.
3. I will acknowledge within 7 days and aim to ship a fix within 30 days of
   triage, depending on severity.

Please give me a chance to ship a fix before public disclosure. Coordinated
disclosure timelines are negotiable for high-severity issues.

## What counts as a vulnerability

- Memory safety in the Cython kernels (buffer over-reads, use-after-free).
- Numerical-correctness bugs that produce silently wrong fits (not the
  documented `rtol~1e-5` drift from R, but cases that violate the M-estimator
  contract).
- Supply-chain issues (e.g., the wheel build pulling untrusted code).

Bugs that produce slow fits, convergence warnings, or fits outside the
documented agreement-with-R tolerances are not security issues; file them as
regular GitHub issues.
