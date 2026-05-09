#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Merge tests/bench/r/*.json + tests/bench/py/*.json into docs/bench-report.md.

For each common case (matched by ``name``), report:

- coefficient agreement: max relative error across all coefficients
- scale agreement: relative error
- covariance agreement: max element-wise relative error on the diagonal
- runtime: R median, py median, ratio (py/R)

Plus a header with versions and platform info.
"""

from __future__ import annotations

import json
import platform
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
R_DIR = REPO_ROOT / "tests" / "bench" / "r"
PY_DIR = REPO_ROOT / "tests" / "bench" / "py"
OUT = REPO_ROOT / "docs" / "bench-report.md"


def _load(path: Path) -> dict:
    return json.loads(path.read_text())


def _max_rel_err(py_coefs: dict, r_coefs: dict) -> float:
    """Max ``|py - R| / max(|R|, 1)`` across all named coefficients."""
    name_map = {"Intercept": "(Intercept)"}
    err = 0.0
    for name, py_v in py_coefs.items():
        r_key = name_map.get(name, name)
        r_v = float(r_coefs.get(r_key, np.nan))
        if not np.isfinite(r_v):
            continue
        err = max(err, abs(py_v - r_v) / max(abs(r_v), 1.0))
    return err


def _diag_rel_err(py_cov: list[list[float]], r_cov: list[list[float]]) -> float:
    py = np.diag(np.asarray(py_cov, dtype=float))
    rr = np.diag(np.asarray(r_cov, dtype=float))
    if py.shape != rr.shape:
        return float("nan")
    denom = np.maximum(np.abs(rr), 1e-12)
    return float(np.max(np.abs(py - rr) / denom))


def _scale_rel_err(py_scale: float, r_scale: float) -> float:
    if r_scale == 0:
        return abs(py_scale - r_scale)
    return abs(py_scale - r_scale) / abs(r_scale)


def _section(title: str, rows: list[list[str]], headers: list[str]) -> str:
    out = [f"## {title}", ""]
    out.append("| " + " | ".join(headers) + " |")
    out.append("|" + "|".join(["---"] * len(headers)) + "|")
    for r in rows:
        out.append("| " + " | ".join(r) + " |")
    out.append("")
    return "\n".join(out)


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)

    r_cases = {p.stem: _load(p) for p in sorted(R_DIR.glob("*.json"))}
    py_cases = {p.stem: _load(p) for p in sorted(PY_DIR.glob("*.json"))}
    common = sorted(set(r_cases) & set(py_cases))
    only_r = sorted(set(r_cases) - set(py_cases))
    only_py = sorted(set(py_cases) - set(r_cases))

    body: list[str] = [
        "# Benchmark report",
        "",
        "Element-wise comparison between `pyrobustlm` and `robustbase::lmrob` "
        "on a fixed corpus of fits. Re-generate with::",
        "",
        "    Rscript scripts/benchmark.R",
        "    python  scripts/benchmark.py",
        "    python  scripts/build_bench_report.py",
        "",
        "## Environment",
        "",
        f"- pyrobustlm: {next(iter(py_cases.values()))['pyrobustlm_version']}",
        f"- Python: {sys.version.split()[0]}",
        f"- Platform: {platform.platform()}",
        f"- robustbase: {next(iter(r_cases.values()))['rb_version']}",
        f"- R: {next(iter(r_cases.values()))['r_version'].split(', ')[0]}",
        "",
    ]

    # ------------------------------------------------------------------
    # Numerical accuracy table
    # ------------------------------------------------------------------
    rows: list[list[str]] = []
    for name in common:
        rj = r_cases[name]
        pj = py_cases[name]
        coef_err = _max_rel_err(pj["coefficients"], rj["coefficients"])
        scale_err = _scale_rel_err(float(pj["scale"]), float(rj["scale"]))
        cov_err = _diag_rel_err(pj["cov"], rj["cov"])
        rows.append(
            [
                name,
                str(rj.get("psi", "")),
                f"{rj['n']}x{rj['p']}",
                f"{coef_err:.2e}",
                f"{scale_err:.2e}",
                f"{cov_err:.2e}",
            ]
        )
    body.append(
        _section(
            "Numerical accuracy: max relative error vs R",
            rows,
            ["case", "psi", "n_x_p", "max coef rerr", "scale rerr", "cov diag max rerr"],
        )
    )

    # ------------------------------------------------------------------
    # Runtime table
    # ------------------------------------------------------------------
    rows = []
    for name in common:
        rj = r_cases[name]
        pj = py_cases[name]
        r_med = float(rj["runtime_median_sec"])
        py_med = float(pj["runtime_median_sec"])
        rows.append(
            [
                name,
                str(rj.get("psi", "")),
                f"{rj['n']}x{rj['p']}",
                f"{1000 * r_med:.1f}",
                f"{1000 * py_med:.1f}",
                f"{py_med / r_med:.2f}x" if r_med > 0 else "n/a",
            ]
        )
    body.append(
        _section(
            "Runtime: median over 5 reps (lower is better)",
            rows,
            ["case", "psi", "n_x_p", "R (ms)", "py (ms)", "py/R"],
        )
    )

    # ------------------------------------------------------------------
    # Coverage / orphans
    # ------------------------------------------------------------------
    body.append("## Coverage")
    body.append("")
    body.append(f"- Cases in both: {len(common)}")
    body.append(f"- Only R: {only_r if only_r else '(none)'}")
    body.append(f"- Only py: {only_py if only_py else '(none)'}")
    body.append("")

    OUT.write_text("\n".join(body))
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
