#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Side-by-side benchmark: default Python+Cython path vs the monolithic
``Control(engine_c=True)`` kernel.

Skips datasets where engine_c lands in a basin that makes vcov_avar1
singular (a known consequence of the Cython subset-draw not being
byte-identical with numpy's choice). Reports median over k_reps runs.
"""

from __future__ import annotations

import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

from pyrobustlm import Control, lmrob

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "tests" / "data"


def _bench(formula: str, df: pd.DataFrame, ctrl: Control, seed: int, k_reps: int = 7) -> float:
    """Return median wall-clock per fit in milliseconds."""
    # Warm up.
    lmrob(formula, df, control=ctrl, seed=seed)
    times = []
    for _ in range(k_reps):
        t0 = time.perf_counter()
        lmrob(formula, df, control=ctrl, seed=seed)
        times.append(time.perf_counter() - t0)
    return float(np.median(times)) * 1000.0


def main() -> None:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

    cases = [
        ("stackloss", "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."),
        ("delivery", "delTime ~ n.prod + distance"),
        ("phosphor", "plant ~ inorg + organic"),
        ("salinity", "Y ~ X1 + X2 + X3"),
    ]
    settings = [(None, "MM"), ("KS2014", "KS2014"), ("KS2011", "KS2011")]

    rows = []
    for dataset, formula in cases:
        path = DATA_DIR / f"{dataset}.csv"
        if not path.exists():
            print(f"skipping {dataset}: data file missing")
            continue
        df = pd.read_csv(path)
        for setting, label in settings:
            try:
                t_def = _bench(formula, df, Control(setting=setting, nResample=500), seed=0)
            except Exception as exc:
                t_def = float("nan")
                print(f"[default] {dataset} {label}: failed ({exc.__class__.__name__})")
            try:
                t_c = _bench(
                    formula,
                    df,
                    Control(setting=setting, nResample=500, engine_c=True),
                    seed=0,
                )
            except Exception as exc:
                t_c = float("nan")
                print(f"[engine_c] {dataset} {label}: failed ({exc.__class__.__name__})")
            speedup = t_def / t_c if t_c and t_c == t_c else float("nan")
            rows.append((dataset, label, t_def, t_c, speedup))

    print()
    print(f"{'dataset':<12} {'setting':<8} {'default(ms)':>13} {'engine_c(ms)':>14} {'speedup':>8}")
    print("-" * 60)
    for dataset, label, t_def, t_c, speedup in rows:
        print(f"{dataset:<12} {label:<8} {t_def:>13.1f} {t_c:>14.1f} {speedup:>7.2f}x")


if __name__ == "__main__":
    main()
