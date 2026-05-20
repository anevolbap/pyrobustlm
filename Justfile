# SPDX-License-Identifier: GPL-3.0-or-later
# Common dev recipes. Install just from https://github.com/casey/just.
# Run `just` (no args) to see the list.

# Show the recipe list.
default:
    @just --list --unsorted

# Full test suite (no coverage gate).
test:
    uv run pytest -q --no-cov

# Full test suite with the coverage gate enforced.
test-cov:
    uv run pytest -q

# Quick filter, e.g. `just t rng`.
t pattern:
    uv run pytest -q --no-cov -k "{{pattern}}"

# Lint + format + type-check.
lint:
    uv run ruff check src tests
    uv run ruff format --check src tests
    uvx ty check src

# Apply ruff format and auto-fixable lint.
fmt:
    uv run ruff format src tests
    uv run ruff check --fix src tests

# Rebuild the Cython extensions in place.
build:
    uv pip install --no-build-isolation -e ".[dev]"

# Build the Sphinx site (fail on warnings, matching RTD).
docs:
    uv run sphinx-build -W -b html docs docs/_build

# Drop the docs build artefacts.
docs-clean:
    rm -rf docs/_build

# Bench against R. Requires Rscript on PATH.
bench:
    Rscript scripts/benchmark.R
    uv run python scripts/benchmark.py
    uv run python scripts/build_bench_report.py

# Release prep: bump _version.py and finalize CHANGELOG. Call as
# `just release 0.5.18`. Doesn't tag; the user reviews + tags by hand.
release version:
    @echo "Bumping to {{version}}"
    sed -i 's/^__version__ = .*/__version__ = "{{version}}"/' src/pylmrob/_version.py
    @echo "Now edit CHANGELOG.md to change '## [Unreleased]' to '## [{{version}}] - $(date +%Y-%m-%d)'"
    @echo "Then: git add -A && git commit -m 'chore(release): v{{version}}' && git tag -a v{{version}} -m 'v{{version}}'"
