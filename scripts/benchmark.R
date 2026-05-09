#!/usr/bin/env Rscript
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Run robustbase::lmrob on a fixed corpus of (dataset, formula, control)
# tuples and dump (coefficients, scale, cov, residuals, weights, runtime)
# to JSON files under tests/bench/r/.
#
# Pair: scripts/benchmark.py runs the same corpus with pyrobustlm and writes
# tests/bench/py/. scripts/build_bench_report.py merges them into
# docs/bench-report.md.

suppressPackageStartupMessages({
  library(robustbase)
  library(jsonlite)
})

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep("--file=", args)
  if (length(match) > 0) {
    return(normalizePath(dirname(sub("--file=", "", args[match[1]]))))
  }
  normalizePath(getwd())
}

REPO_ROOT <- normalizePath(file.path(script_dir(), ".."))
OUT_DIR <- file.path(REPO_ROOT, "tests", "bench", "r")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Corpus: classical + synthetic + per-psi-family
# ---------------------------------------------------------------------------
classical <- list(
  list(name="classical_stackloss", dataset="stackloss",
       formula = stack.loss ~ Air.Flow + Water.Temp + Acid.Conc., control=list()),
  list(name="classical_coleman", dataset="coleman",  formula = Y ~ ., control=list()),
  list(name="classical_salinity", dataset="salinity", formula = Y ~ ., control=list()),
  list(name="classical_delivery", dataset="delivery", formula = delTime ~ n.prod + distance, control=list()),
  list(name="classical_phosphor", dataset="phosphor", formula = plant ~ inorg + organic, control=list()),
  list(name="classical_aircraft", dataset="aircraft", formula = Y ~ X1 + X2 + X3 + X4, control=list()),
  list(name="classical_pension",  dataset="pension",  formula = Reserves ~ Income, control=list()),
  list(name="classical_starsCYG", dataset="starsCYG", formula = log.light ~ log.Te, control=list()),
  list(name="classical_hbk",      dataset="hbk",      formula = Y ~ ., control=list()),
  list(name="classical_wood",     dataset="wood",     formula = y ~ ., control=list())
)

per_psi <- list(
  list(name="setting_KS2014_stackloss", dataset="stackloss",
       formula=stack.loss ~ Air.Flow + Water.Temp + Acid.Conc., control=list(setting="KS2014")),
  list(name="setting_KS2011_stackloss", dataset="stackloss",
       formula=stack.loss ~ Air.Flow + Water.Temp + Acid.Conc., control=list(setting="KS2011")),
  list(name="psi_bisquare", dataset="stackloss",
       formula=stack.loss ~ Air.Flow + Water.Temp + Acid.Conc., control=list(psi="bisquare")),
  list(name="psi_optimal",  dataset="stackloss",
       formula=stack.loss ~ Air.Flow + Water.Temp + Acid.Conc., control=list(psi="optimal")),
  list(name="psi_hampel",   dataset="stackloss",
       formula=stack.loss ~ Air.Flow + Water.Temp + Acid.Conc., control=list(psi="hampel")),
  list(name="psi_lqq",      dataset="stackloss",
       formula=stack.loss ~ Air.Flow + Water.Temp + Acid.Conc., control=list(psi="lqq")),
  list(name="psi_ggw",      dataset="stackloss",
       formula=stack.loss ~ Air.Flow + Water.Temp + Acid.Conc., control=list(psi="ggw"))
)

# Synthetic timing grid: focus on speed comparison
make_synthetic <- function(n, p, contam=0.1, seed=1) {
  set.seed(seed)
  X <- matrix(rnorm(n*p), n, p)
  beta <- runif(p, -2, 2)
  y <- as.numeric(X %*% beta) + rnorm(n)
  if (contam > 0) {
    n_bad <- ceiling(n * contam)
    idx <- sample.int(n, n_bad)
    y[idx] <- y[idx] + rnorm(n_bad, 20, 1)
  }
  df <- as.data.frame(X)
  names(df) <- paste0("x", seq_len(p))
  df$y <- y
  df
}

synthetic <- list(
  list(name="synth_n100_p5",     n=100,   p=5),
  list(name="synth_n500_p10",    n=500,   p=10),
  list(name="synth_n1000_p10",   n=1000,  p=10),
  list(name="synth_n2000_p20",   n=2000,  p=20),
  list(name="synth_n5000_p20",   n=5000,  p=20),
  list(name="synth_n10000_p20",  n=10000, p=20),
  list(name="synth_n10000_p50",  n=10000, p=50)
)

# Per-psi-family timings for a moderate-sized synthetic, so we can
# compare runtime profile across families on more than just stackloss.
synth_per_psi <- list()
for (fam in c("bisquare", "optimal", "hampel", "lqq", "ggw")) {
  for (cfg in list(c(500, 10), c(2000, 20))) {
    n <- cfg[1]; p <- cfg[2]
    synth_per_psi[[length(synth_per_psi) + 1]] <- list(
      name = sprintf("synth_%s_n%d_p%d", fam, n, p),
      n = n, p = p, psi_family = fam
    )
  }
}

corpus <- c(classical, per_psi)

# ---------------------------------------------------------------------------
load_dataset <- function(entry) {
  if (!is.null(entry$synth_args)) {
    df <- do.call(make_synthetic, entry$synth_args)
    f <- as.formula(paste("y ~", paste(setdiff(names(df), "y"), collapse=" + ")))
    return(list(data=df, formula=f))
  }
  e <- new.env()
  for (pkg in c("robustbase", "datasets")) {
    ok <- suppressWarnings(tryCatch({
      data(list = entry$dataset, package = pkg, envir = e)
      exists(entry$dataset, envir = e)
    }, error = function(...) FALSE))
    if (isTRUE(ok)) break
  }
  list(data = get(entry$dataset, envir = e), formula = entry$formula)
}

run_one <- function(entry) {
  cat(sprintf("[R bench] %s\n", entry$name))
  ds <- load_dataset(entry)
  ctrl <- do.call(lmrob.control, c(list(nResample = 500), entry$control))
  set.seed(1)
  # Warm-up
  fit_w <- tryCatch(lmrob(ds$formula, data=ds$data, control=ctrl), error=function(e) NULL)
  if (is.null(fit_w)) {
    cat("  ERROR\n"); return(invisible(NULL))
  }
  # Time over k repetitions
  k <- 5
  set.seed(1)
  t <- numeric(k)
  for (i in 1:k) {
    set.seed(1)
    tt <- system.time(fit <- lmrob(ds$formula, data=ds$data, control=ctrl))
    t[i] <- as.numeric(tt["elapsed"])
  }
  out <- list(
    name = entry$name,
    coefficients = as.list(coef(fit)),
    scale = fit$scale,
    cov = unname(as.matrix(vcov(fit))),
    converged = isTRUE(fit$converged),
    rweights = unname(fit$rweights),
    residuals = unname(residuals(fit)),
    fitted = unname(fitted(fit)),
    psi = fit$control$psi,
    n = nrow(ds$data),
    p = length(coef(fit)),
    runtimes_sec = t,
    runtime_min_sec = min(t),
    runtime_median_sec = median(t),
    rb_version = as.character(packageVersion("robustbase")),
    r_version = R.version.string
  )
  write(jsonlite::toJSON(out, auto_unbox=TRUE, digits=17, na="null", null="null"),
        file = file.path(OUT_DIR, paste0(entry$name, ".json")))
}

run_synthetic <- function(s, psi_family = NULL) {
  ctrl <- if (is.null(psi_family)) list() else list(psi = psi_family)
  entry <- list(name=s$name, dataset=NULL, formula=NULL, control=ctrl,
                synth_args=list(n=s$n, p=s$p))
  run_one(entry)
}

for (e in corpus) run_one(e)
for (s in synthetic) run_synthetic(s)
for (s in synth_per_psi) run_synthetic(s, psi_family = s$psi_family)
cat(sprintf("Wrote R bench JSONs to %s\n", OUT_DIR))
