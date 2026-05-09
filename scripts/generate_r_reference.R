#!/usr/bin/env Rscript
# SPDX-License-Identifier: GPL-3.0-or-later
#
# scripts/generate_r_reference.R
#
# Run R's robustbase::lmrob on a fixed corpus of (dataset, formula, control,
# seed) tuples and dump the results to JSON files under tests/reference/.
#
# Usage:
#   Rscript scripts/generate_r_reference.R
#
# Re-run whenever the corpus changes. JSON files are committed so CI can
# diff Python output against R offline.

suppressPackageStartupMessages({
  library(robustbase)
  library(jsonlite)
})

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  match <- grep(file_arg, args)
  if (length(match) > 0) {
    return(normalizePath(dirname(sub(file_arg, "", args[match[1]]))))
  }
  if (interactive() && !is.null(sys.frame(1)$ofile)) {
    return(normalizePath(dirname(sys.frame(1)$ofile)))
  }
  normalizePath(getwd())
}

REPO_ROOT <- normalizePath(file.path(script_dir(), ".."))
OUT_DIR   <- file.path(REPO_ROOT, "tests", "reference")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Corpus
# ---------------------------------------------------------------------------
# Each entry: list(name, dataset, formula, control_args, seed).
# Datasets are loaded from robustbase via data().
corpus <- list(
  list(
    name    = "stackloss_default",
    dataset = "stackloss",
    formula = stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.,
    control = list(),
    seed    = 1
  ),
  list(
    name    = "coleman_default",
    dataset = "coleman",
    formula = Y ~ .,
    control = list(),
    seed    = 1
  ),
  list(
    name    = "salinity_default",
    dataset = "salinity",
    formula = Y ~ .,
    control = list(),
    seed    = 1
  ),
  list(
    name    = "wood_default",
    dataset = "wood",
    formula = y ~ .,
    control = list(),
    seed    = 1
  ),
  list(
    name    = "hbk_default",
    dataset = "hbk",
    formula = Y ~ .,
    control = list(),
    seed    = 1
  ),
  list(
    name    = "starsCYG_default",
    dataset = "starsCYG",
    formula = log.light ~ log.Te,
    control = list(),
    seed    = 1
  ),
  list(
    name    = "stackloss_KS2011",
    dataset = "stackloss",
    formula = stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.,
    control = list(setting = "KS2011"),
    seed    = 1
  ),
  list(
    name    = "stackloss_KS2014",
    dataset = "stackloss",
    formula = stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.,
    control = list(setting = "KS2014"),
    seed    = 1
  ),
  list(
    name    = "stackloss_psi_bisquare",
    dataset = "stackloss",
    formula = stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.,
    control = list(psi = "bisquare"),
    seed    = 1
  ),
  list(
    name    = "stackloss_psi_optimal",
    dataset = "stackloss",
    formula = stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.,
    control = list(psi = "optimal"),
    seed    = 1
  ),
  list(
    name    = "stackloss_psi_ggw",
    dataset = "stackloss",
    formula = stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.,
    control = list(psi = "ggw"),
    seed    = 1
  ),
  list(
    name    = "stackloss_psi_lqq",
    dataset = "stackloss",
    formula = stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.,
    control = list(psi = "lqq"),
    seed    = 1
  ),
  # NOTE: huber not added because robustbase::lmrob requires a redescending
  # psi for the MM step. huber stays available via psi="huber" in our M and
  # M-scale APIs.
  list(
    name    = "stackloss_psi_hampel",
    dataset = "stackloss",
    formula = stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.,
    control = list(psi = "hampel"),
    seed    = 1
  ),
  list(
    name    = "delivery_default",
    dataset = "delivery",
    formula = delTime ~ n.prod + distance,
    control = list(),
    seed    = 1
  ),
  list(
    name    = "aircraft_default",
    dataset = "aircraft",
    formula = Y ~ X1 + X2 + X3 + X4,
    control = list(),
    seed    = 1
  ),
  list(
    name    = "pension_default",
    dataset = "pension",
    formula = Reserves ~ Income,
    control = list(),
    seed    = 1
  ),
  list(
    name    = "phosphor_default",
    dataset = "phosphor",
    formula = plant ~ inorg + organic,
    control = list(),
    seed    = 1
  ),
  # education has Region (factor) and X1..X3 (numeric). State (50 levels) is
  # left out for now; Phase 5 (M-S init) will add a full-State reference.
  list(
    name    = "education_simple",
    dataset = "education",
    formula = Y ~ Region + X1 + X2 + X3,
    control = list(),
    seed    = 1
  )
)

# ---------------------------------------------------------------------------
# Synthetic configurations
# ---------------------------------------------------------------------------
make_synthetic <- function(n, p, contam_frac, seed = 1) {
  set.seed(seed)
  X    <- matrix(rnorm(n * p), n, p)
  beta <- runif(p, -2, 2)
  eps  <- rnorm(n)
  y    <- as.numeric(X %*% beta) + eps
  if (contam_frac > 0) {
    n_bad <- ceiling(n * contam_frac)
    idx   <- sample.int(n, n_bad)
    y[idx] <- y[idx] + rnorm(n_bad, mean = 20, sd = 1)
  }
  df         <- as.data.frame(X)
  names(df)  <- paste0("x", seq_len(p))
  df$y       <- y
  attr(df, "true_beta") <- beta
  df
}

synth_grid <- expand.grid(
  n           = c(50, 200, 1000),
  p           = c(3, 10),
  contam_frac = c(0.0, 0.10, 0.30)
)
for (i in seq_len(nrow(synth_grid))) {
  row <- synth_grid[i, ]
  name <- sprintf("synthetic_n%d_p%d_c%02d",
                  row$n, row$p, as.integer(row$contam_frac * 100))
  corpus[[length(corpus) + 1]] <- list(
    name    = name,
    dataset = paste0("__synthetic_", i),
    formula = NULL,
    control = list(),
    seed    = 1L,
    synth   = list(n = row$n, p = row$p, contam_frac = row$contam_frac)
  )
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
load_dataset <- function(entry) {
  if (!is.null(entry$synth)) {
    df <- make_synthetic(entry$synth$n, entry$synth$p,
                         entry$synth$contam_frac, entry$seed)
    f  <- as.formula(paste("y ~", paste(setdiff(names(df), "y"), collapse = " + ")))
    return(list(data = df, formula = f))
  }
  e <- new.env()
  # Try robustbase first, then base 'datasets' (e.g. stackloss).
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
  cat(sprintf("[ref] %s\n", entry$name))
  ds <- load_dataset(entry)
  ctrl <- do.call(lmrob.control, entry$control)
  set.seed(entry$seed)
  fit <- tryCatch(
    lmrob(ds$formula, data = ds$data, control = ctrl),
    error = function(e) {
      cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(fit)) return(invisible(NULL))

  out <- list(
    name        = entry$name,
    dataset     = entry$dataset,
    formula     = deparse(ds$formula),
    control     = entry$control,
    seed        = entry$seed,
    coefficients = as.list(coef(fit)),
    scale        = fit$scale,
    weights      = unname(fit$rweights),
    residuals    = unname(residuals(fit)),
    fitted       = unname(fitted(fit)),
    cov          = unname(as.matrix(vcov(fit))),
    df_residual  = fit$df.residual,
    converged    = isTRUE(fit$converged),
    init_S       = list(
      coefficients = if (!is.null(fit$init.S)) as.list(fit$init.S$coefficients) else NULL,
      scale        = if (!is.null(fit$init.S)) fit$init.S$scale else NULL
    ),
    rweights     = unname(fit$rweights),
    psi          = fit$control$psi,
    tuning_psi   = fit$control$tuning.psi,
    tuning_chi   = fit$control$tuning.chi,
    rb_version   = as.character(packageVersion("robustbase")),
    r_version    = R.version.string
  )

  out_path <- file.path(OUT_DIR, paste0(entry$name, ".json"))
  write(jsonlite::toJSON(out, auto_unbox = TRUE, digits = 17,
                         na = "null", null = "null"),
        file = out_path)
}

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
for (entry in corpus) run_one(entry)
cat(sprintf("Wrote references to %s\n", OUT_DIR))
