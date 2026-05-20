# Summary

<!-- One paragraph describing what changes and why. -->

# Test plan

- [ ] Full local suite: `pytest --no-cov` passes.
- [ ] Ruff + ty clean: `ruff check src tests && ruff format --check src tests && uvx ty check src`.
- [ ] Touched `CHANGELOG.md` under `[Unreleased]`, or this PR is doc-only / pure chore.
- [ ] If touching `src/pylmrob/_core/*.pyx`: rebuilt via `uv pip install --no-build-isolation -e ".[dev]"` and reran the suite.

# Notes

<!-- Anything the reviewer needs to know that isn't obvious from the diff: -->
<!-- - Numerical agreement against R, if applicable. -->
<!-- - Performance impact, if applicable. -->
<!-- - Cross-references to issues, PRs, or external sources. -->
