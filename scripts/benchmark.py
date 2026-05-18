#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Run pylmrob on the same corpus as scripts/benchmark.R.

Writes ``tests/bench/py/<name>.json`` for each fit with coefficients,
scale, covariance, residuals, weights, and timing data. Pair with
``scripts/build_bench_report.py`` to produce ``docs/bench-report.md``.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

from pylmrob import Control, lmrob

REPO_ROOT = Path(__file__).resolve().parents[1]

_ENGINE_C = False
OUT_DIR = REPO_ROOT / "tests" / "bench" / "py"


def _ensure_dataset(name: str) -> pd.DataFrame:
    out = REPO_ROOT / "tests" / "data" / f"{name}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    if not out.exists():
        subprocess.run(
            [
                "Rscript",
                "-e",
                f"library(robustbase); data({name}); write.csv({name}, '{out}', row.names=FALSE)",
            ],
            capture_output=True,
            check=True,
        )
    return pd.read_csv(out)


def _make_synthetic(n: int, p: int, contam: float = 0.1, seed: int = 1) -> pd.DataFrame:
    """Mirror the R-side synthetic generator (uses R's RNG for parity)."""
    out = REPO_ROOT / "tests" / "data" / f"synth_n{n}_p{p}.csv"
    if not out.exists():
        # Generate via R so the Python and R fits see the *same* data.
        rscript = (
            f"set.seed(1);"
            f"X <- matrix(rnorm({n}*{p}), {n}, {p}); "
            f"beta <- runif({p}, -2, 2); "
            f"y <- as.numeric(X %*% beta) + rnorm({n}); "
            f"n_bad <- ceiling({n}*{contam}); "
            f"idx <- sample.int({n}, n_bad); "
            f"y[idx] <- y[idx] + rnorm(n_bad, 20, 1); "
            f"df <- as.data.frame(X); "
            f"names(df) <- paste0('x', seq_len({p})); "
            f"df$y <- y; "
            f"write.csv(df, '{out}', row.names=FALSE)"
        )
        subprocess.run(["Rscript", "-e", rscript], capture_output=True, check=True)
    return pd.read_csv(out)


CLASSICAL = [
    (
        "classical_stackloss",
        "stackloss",
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        "bisquare",
    ),
    ("classical_coleman", "coleman", "Y ~ .", "bisquare"),
    ("classical_salinity", "salinity", "Y ~ .", "bisquare"),
    ("classical_delivery", "delivery", "delTime ~ n.prod + distance", "bisquare"),
    ("classical_phosphor", "phosphor", "plant ~ inorg + organic", "bisquare"),
    ("classical_aircraft", "aircraft", "Y ~ X1 + X2 + X3 + X4", "bisquare"),
    ("classical_pension", "pension", "Reserves ~ Income", "bisquare"),
    ("classical_starsCYG", "starsCYG", "log.light ~ log.Te", "bisquare"),
    ("classical_hbk", "hbk", "Y ~ .", "bisquare"),
    ("classical_wood", "wood", "y ~ .", "bisquare"),
]

PER_PSI = [
    ("psi_bisquare", "bisquare"),
    ("psi_optimal", "optimal"),
    ("psi_hampel", "hampel"),
    ("psi_lqq", "lqq"),
    ("psi_ggw", "ggw"),
]

PER_SETTING = [
    ("setting_KS2014_stackloss", "KS2014"),
    ("setting_KS2011_stackloss", "KS2011"),
]

SYNTHETIC = [
    ("synth_n100_p5", 100, 5),
    ("synth_n500_p10", 500, 10),
    ("synth_n1000_p10", 1000, 10),
    ("synth_n2000_p20", 2000, 20),
    ("synth_n5000_p20", 5000, 20),
    ("synth_n10000_p20", 10000, 20),
    ("synth_n10000_p50", 10000, 50),
]

SYNTH_PER_PSI = []
for fam in ("bisquare", "optimal", "hampel", "lqq", "ggw"):
    for n, p in [(500, 10), (2000, 20)]:
        SYNTH_PER_PSI.append((f"synth_{fam}_n{n}_p{p}", n, p, fam))


def _fit_and_time(
    formula: str,
    df: pd.DataFrame,
    psi_family: str = "bisquare",
    setting: str | None = None,
    k_reps: int = 5,
) -> dict:
    # ``setting=None`` matches R's plain ``lmrob.control()`` default
    # (psi=bisquare unless overridden, method=MM, cov=.vcov.avar1).
    ctrl = Control(psi=psi_family, setting=setting, nResample=500, engine_c=_ENGINE_C)
    _ = lmrob(formula, df, control=ctrl, seed=1)
    runtimes = []
    for _ in range(k_reps):
        t0 = time.perf_counter()
        fit = lmrob(formula, df, control=ctrl, seed=1)
        runtimes.append(time.perf_counter() - t0)
    return {
        "fit": fit,
        "runtimes_sec": runtimes,
        "runtime_min_sec": float(min(runtimes)),
        "runtime_median_sec": float(np.median(runtimes)),
    }


def _serialize(result: dict, name: str, psi_family: str) -> None:
    fit = result["fit"]
    name_map = {"Intercept": "(Intercept)"}
    coefs = {name_map.get(n, n): float(v) for n, v in zip(fit.term_names_, fit.coef_, strict=True)}
    out = {
        "name": name,
        "coefficients": coefs,
        "scale": float(fit.scale_),
        "cov": fit.cov_.tolist(),
        "converged": bool(fit.converged_),
        "rweights": fit.rweights_.tolist(),
        "residuals": fit.residuals_.tolist(),
        "fitted": fit.fitted_.tolist(),
        "psi": psi_family,
        "n": int(fit.nobs_),
        "p": int(fit.coef_.size),
        "runtimes_sec": result["runtimes_sec"],
        "runtime_min_sec": result["runtime_min_sec"],
        "runtime_median_sec": result["runtime_median_sec"],
        "py_version": __import__("sys").version.split()[0],
        "pylmrob_version": __import__("pylmrob").__version__,
    }
    (OUT_DIR / f"{name}.json").write_text(json.dumps(out, indent=2))


def _safe_fit(name: str, fn, *args, **kwargs):
    try:
        return fn(*args, **kwargs)
    except Exception as exc:
        print(f"  SKIP {name}: {exc.__class__.__name__}: {exc}")
        return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--engine-c",
        action="store_true",
        help="Run with Control(engine_c=True); write to tests/bench/py_engine_c/.",
    )
    args = parser.parse_args()
    global _ENGINE_C, OUT_DIR
    _ENGINE_C = args.engine_c
    OUT_DIR = REPO_ROOT / "tests" / "bench" / ("py_engine_c" if _ENGINE_C else "py")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"[py bench] engine_c={_ENGINE_C} writing to {OUT_DIR.relative_to(REPO_ROOT)}")

    for name, dataset, formula, psi_family in CLASSICAL:
        print(f"[py bench] {name}")
        df = _ensure_dataset(dataset)
        if formula == "Y ~ .":
            rhs = " + ".join([c for c in df.columns if c != "Y"])
            formula = f"Y ~ {rhs}"
        elif formula == "y ~ .":
            rhs = " + ".join([c for c in df.columns if c != "y"])
            formula = f"y ~ {rhs}"
        result = _safe_fit(name, _fit_and_time, formula, df, psi_family)
        if result is not None:
            _serialize(result, name, psi_family)

    for name, psi_family in PER_PSI:
        print(f"[py bench] {name}")
        df = _ensure_dataset("stackloss")
        result = _safe_fit(
            name,
            _fit_and_time,
            "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
            df,
            psi_family,
        )
        if result is not None:
            _serialize(result, name, psi_family)

    def _per_setting(name, setting):
        df = _ensure_dataset("stackloss")
        ctrl = Control(setting=setting, nResample=500, engine_c=_ENGINE_C)
        _ = lmrob(
            "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
            df,
            control=ctrl,
            seed=1,
        )
        runtimes = []
        for _ in range(5):
            t0 = time.perf_counter()
            fit = lmrob(
                "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
                df,
                control=ctrl,
                seed=1,
            )
            runtimes.append(time.perf_counter() - t0)
        return {
            "fit": fit,
            "runtimes_sec": runtimes,
            "runtime_min_sec": float(min(runtimes)),
            "runtime_median_sec": float(np.median(runtimes)),
        }, ctrl.psi

    for name, setting in PER_SETTING:
        print(f"[py bench] {name}")
        out = _safe_fit(name, _per_setting, name, setting)
        if out is not None:
            result, psi = out
            _serialize(result, name, psi)

    for name, n, p in SYNTHETIC:
        print(f"[py bench] {name}")
        df = _make_synthetic(n, p)
        formula = "y ~ " + " + ".join(f"x{i + 1}" for i in range(p))
        result = _safe_fit(name, _fit_and_time, formula, df, "bisquare")
        if result is not None:
            _serialize(result, name, "bisquare")

    for name, n, p, fam in SYNTH_PER_PSI:
        print(f"[py bench] {name}")
        df = _make_synthetic(n, p)
        formula = "y ~ " + " + ".join(f"x{i + 1}" for i in range(p))
        result = _safe_fit(name, _fit_and_time, formula, df, fam)
        if result is not None:
            _serialize(result, name, fam)
    print(f"Wrote Python bench JSONs to {OUT_DIR}")


if __name__ == "__main__":
    main()
