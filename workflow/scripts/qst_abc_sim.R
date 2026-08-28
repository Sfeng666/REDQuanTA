#!/usr/bin/env Rscript
# ABC-based QST estimation - adapted from run_abc_sim.R
# MEMORY-OPTIMIZED VERSION: Processes in batches to reduce peak memory usage
#
# This script estimates QST using ABC with 100,000 simulations
# Can be used for both trait QST and neutral QST estimation
#
# Usage:
#   Rscript qst_abc_sim.R trait obs_stats.RData ignored output.RData [num_sim] [summary_stats]
#   Rscript qst_abc_sim.R neutral 0.05 0.1 output.RData [num_sim] [summary_stats]
#
# Arguments:
#   1. mode: "trait" or "neutral"
#   2. input_file: For trait mode: path to observed summary stats
#                  For neutral mode: FST value
#   3. ext_sd: Extrinsic SD (from observed trait data)
#   4. output_file: Path to save results
#   5. num_sim: Number of ABC simulations (default: 100000)
#   6. summary_stats: Comma-separated summary statistics (default: "QST,ratioVbetweenVtotal")

suppressPackageStartupMessages({
  library(abc)
  library(parallel)
})

# Fast loclinear / rejection (no abc::abc for those methods); same directory as this script.
{
  e <- Sys.getenv("QST_ABC_SCRIPT_DIR", unset = "")
  script_dir <- if (nzchar(e)) {
    e
  } else {
    a <- commandArgs(trailingOnly = FALSE)
    f <- grep("^--file=", a, value = TRUE)
    if (length(f)) {
      dirname(normalizePath(sub("^--file=", "", f[[1]]), winslash = "/"))
    } else {
      getwd()
    }
  }
  fast_src <- file.path(script_dir, "abc_fast_estimators.R")
  if (!file.exists(fast_src)) stop("Missing required file: ", fast_src)
  sys.source(fast_src, envir = environment())
}

# QST_ABC_METHOD: neuralnet | rejection | loclinear | ridge (rej -> rejection). Default loclinear.
# loclinear and rejection use the vectorized implementation in abc_fast_estimators.R (same
# numerical results as abc::abc on a fixed reference pool). neuralnet and ridge use abc::abc.
get_abc_method <- function() {
  m <- trimws(Sys.getenv("QST_ABC_METHOD", "loclinear"))
  if (m == "rej") m <- "rejection"
  if (!nzchar(m)) m <- "loclinear"
  m
}

# Relative ABC tolerance: tol = max(0.001, QST_ABC_TOL_NUMERATOR / n_valid). Default 50 (matches legacy).
get_abc_tol_numerator <- function() {
  v <- suppressWarnings(as.integer(Sys.getenv("QST_ABC_TOL_NUMERATOR", "50")))
  if (length(v) != 1L || is.na(v) || v < 1L) 50L else v
}
run_abc_qst <- function(target, param, sumstat, tol, transf) {
  meth <- get_abc_method()
  if (meth %in% c("loclinear", "rejection")) {
    res <- fast_abc_estimate_one(
      target, param, sumstat, tol,
      estimator = meth, hcorr = TRUE
    )
    return(list(
      fast_scalar_qst = res$qst,
      abc_among_pop_sd = res$abc_among_pop_sd,
      abc_within_pop_sd = res$abc_within_pop_sd,
      abc_ext_sd = res$abc_ext_sd,
      qst_lo = res$qst_lo,
      qst_hi = res$qst_hi,
      ratioVext_lo = res$ratioVext_lo,
      ratioVext_hi = res$ratioVext_hi
    ))
  }
  if (meth == "neuralnet") {
    res <- abc(
      target = target, param = param, sumstat = sumstat, tol = tol, transf = transf,
      method = meth, sizenet = 10
    )
    return(abc_result_with_quantiles(res))
  }
  abc_result_with_quantiles(abc(
    target = target, param = param, sumstat = sumstat, tol = tol, transf = transf, method = meth
  ))
}

abc_result_with_quantiles <- function(res_abc) {
  m <- res_abc$adj.values
  if (is.null(m) || !is.matrix(m) || nrow(m) < 1L) m <- res_abc$unadj.values
  if (is.null(m) || !is.matrix(m) || nrow(m) < 1L) {
    return(list(
      fast_scalar_qst = NA_real_,
      abc_among_pop_sd = NA_real_,
      abc_within_pop_sd = NA_real_,
      abc_ext_sd = NA_real_,
      qst_lo = NA_real_,
      qst_hi = NA_real_,
      ratioVext_lo = NA_real_,
      ratioVext_hi = NA_real_
    ))
  }
  q <- abc_posterior_quantiles_from_param(m)
  sds <- abc_sd_means_from_param(m)
  c(
    list(fast_scalar_qst = qst_from_param_matrix(m)),
    sds,
    list(
      qst_lo = q$qst_lo,
      qst_hi = q$qst_hi,
      ratioVext_lo = q$ratioVext_lo,
      ratioVext_hi = q$ratioVext_hi
    )
  )
}
qst_mean_from_abc <- function(res_abc) {
  if (is.list(res_abc) && !is.null(res_abc$fast_scalar_qst)) {
    v <- res_abc$fast_scalar_qst
    if (length(v) == 1L && is.finite(v)) return(as.numeric(v))
    return(NA_real_)
  }
  m <- res_abc$adj.values
  if (is.null(m) || !is.matrix(m) || nrow(m) < 1L) m <- res_abc$unadj.values
  if (is.null(m) || !is.matrix(m) || nrow(m) < 1L) return(NA_real_)
  mean(m[, "sd_between_pop"]^2 / (m[, "sd_between_pop"]^2 + 2 * m[, "sd_within_pop"]^2))
}

abc_sd_means_from_abc <- function(res_abc) {
  if (is.list(res_abc) && !is.null(res_abc$abc_among_pop_sd)) {
    return(list(
      abc_among_pop_sd = as.numeric(res_abc$abc_among_pop_sd),
      abc_within_pop_sd = as.numeric(res_abc$abc_within_pop_sd),
      abc_ext_sd = as.numeric(res_abc$abc_ext_sd)
    ))
  }
  m <- res_abc$adj.values
  if (is.null(m) || !is.matrix(m) || nrow(m) < 1L) m <- res_abc$unadj.values
  if (is.null(m) || !is.matrix(m) || nrow(m) < 1L) {
    return(list(
      abc_among_pop_sd = NA_real_,
      abc_within_pop_sd = NA_real_,
      abc_ext_sd = NA_real_
    ))
  }
  abc_sd_means_from_param(m)
}

abc_na_result <- function() {
  list(
    qst = NA_real_,
    qst_lo = NA_real_,
    qst_hi = NA_real_,
    abc_among_pop_sd = NA_real_,
    abc_within_pop_sd = NA_real_,
    abc_ext_sd = NA_real_,
    abc_ratioVext = NA_real_,
    abc_ratioVext_lo = NA_real_,
    abc_ratioVext_hi = NA_real_
  )
}

abc_sds_to_ratioVext <- function(among_sd, within_sd, ext_sd) {
  if (!all(is.finite(c(among_sd, within_sd, ext_sd)))) return(NA_real_)
  vg <- (among_sd^2 + 2 * within_sd^2) / 2
  if (!is.finite(vg) || vg <= 0) return(NA_real_)
  ext_sd^2 / vg
}

ratioVext_from_trait_qst_rdata <- function(path) {
  if (!file.exists(path)) stop("Missing trait QST file: ", path)
  env <- new.env()
  load(path, envir = env)
  if (!exists("result", envir = env)) stop("No result object in ", path)
  r <- env$result
  if (isTRUE(as.integer(r$both_var_negative) == 1L)) {
    return(list(ratioVext = NA_real_, source = "abc_trait", valid = FALSE))
  }
  rv <- abc_sds_to_ratioVext(r$abc_among_pop_sd, r$abc_within_pop_sd, r$abc_ext_sd)
  list(ratioVext = rv, source = "abc_trait", valid = is.finite(rv))
}

resolve_ratioVext_arg <- function(ratioVext_arg) {
  if (file.exists(ratioVext_arg) && grepl("_trait_qst\\.RData$", ratioVext_arg, ignore.case = TRUE)) {
    info <- ratioVext_from_trait_qst_rdata(ratioVext_arg)
    cat("RatioVext from ABC trait_qst:", info$ratioVext, "(source:", info$source, ")\n")
    return(info)
  }
  if (file.exists(ratioVext_arg)) {
    load(ratioVext_arg)
    rv <- as.numeric(obs_stats["ratioVext"])
    cat("RatioVext from obs_stats file:", rv, "\n")
    return(list(ratioVext = rv, source = "anova_obs_stats", valid = is.finite(rv)))
  }
  rv <- as.numeric(ratioVext_arg)
  cat("RatioVext from argument:", rv, "\n")
  list(ratioVext = rv, source = "numeric", valid = is.finite(rv))
}

abc_result_to_estimate <- function(res_abc) {
  sds <- abc_sd_means_from_abc(res_abc)
  qst <- qst_mean_from_abc(res_abc)
  qlo <- if (!is.null(res_abc$qst_lo)) res_abc$qst_lo else NA_real_
  qhi <- if (!is.null(res_abc$qst_hi)) res_abc$qst_hi else NA_real_
  rvlo <- if (!is.null(res_abc$ratioVext_lo)) res_abc$ratioVext_lo else NA_real_
  rvhi <- if (!is.null(res_abc$ratioVext_hi)) res_abc$ratioVext_hi else NA_real_
  rv <- abc_sds_to_ratioVext(sds$abc_among_pop_sd, sds$abc_within_pop_sd, sds$abc_ext_sd)
  list(
    qst = qst,
    qst_lo = qlo,
    qst_hi = qhi,
    abc_among_pop_sd = sds$abc_among_pop_sd,
    abc_within_pop_sd = sds$abc_within_pop_sd,
    abc_ext_sd = sds$abc_ext_sd,
    abc_ratioVext = rv,
    abc_ratioVext_lo = rvlo,
    abc_ratioVext_hi = rvhi
  )
}

# Load e1071 only if needed (for skewness and kurtosis)
.e1071_loaded <- FALSE
load_e1071_if_needed <- function() {
  if (!.e1071_loaded) {
    suppressPackageStartupMessages(library(e1071))
    .e1071_loaded <<- TRUE
  }
}

# Use fewer cores to limit memory usage from forking
# Each forked process duplicates memory, so fewer cores = lower peak memory
# 4 cores is a good balance between speed and memory for CHTC jobs with 4GB limit
num_cores <- min(4, max(1, detectCores() - 1))

# All possible summary statistics (expanded pool; skewness/kurtosis removed from default pool)
ALL_SUMMARY_STATS <- c(
  "among_pop_sd", "within_pop_sd", "ext_sd", "QST", "ratioVext",
  "F_among_pop", "F_within_pop",
  "ratioVbetweenVext", "ratioVbetweenVtotal", "ratioVwithinVtotal", "ratioVextVtotal",
  "both_var_negative"
)

# Basic stats always calculated (needed for variance components)
BASIC_STATS <- c("among_pop_sd", "within_pop_sd", "ext_sd", "both_var_negative")

#' True when both ANOVA genetic variance components are non-positive (QST not estimable).
#' Uses only the explicit both_var_negative flag from the floor policy.
#' (Do not treat among_pop_sd == within_pop_sd as a proxy; that is wrong under the ridge floor.)
is_both_neg <- function(obs_stats) {
  if ("both_var_negative" %in% names(obs_stats)) {
    v <- suppressWarnings(as.numeric(obs_stats["both_var_negative"]))
    if (length(v) == 1L && !is.na(v) && v == 1) return(TRUE)
  }
  FALSE
}

#' Determine which summary statistics need to be calculated
#' Based on the combinations that will be used for ABC
#' 
#' @param summary_stat_combos List of character vectors specifying combinations
#' @return Character vector of unique stats that need to be calculated
get_required_stats <- function(summary_stat_combos) {
  # Get all unique stats from all combinations
  all_stats <- unique(unlist(summary_stat_combos))
  
  # Always include basic stats (needed for prior generation and variance calculation)
  required <- union(BASIC_STATS, all_stats)
  
  # Order according to ALL_SUMMARY_STATS for consistent column ordering
  required <- intersect(ALL_SUMMARY_STATS, required)
  
  return(required)
}

#' Parse summary stats argument
#' Supports: comma-separated string or file path
#' 
#' @param summary_stats_arg The argument value (string, "all", or file path)
#' @return List of character vectors, each specifying a combination
parse_summary_stats <- function(summary_stats_arg) {
  # Check if it's a file path
  if (file.exists(summary_stats_arg)) {
    cat("Loading summary stat combinations from file:", summary_stats_arg, "\n")
    lines <- readLines(summary_stats_arg)
    lines <- lines[lines != ""]  # Remove empty lines
    
    # Parse tab-separated combinations
    combos <- lapply(lines, function(line) {
      stats <- strsplit(line, "\t")[[1]]
      # Handle ratioVenv -> ratioVext conversion
      stats <- gsub("ratioVenv", "ratioVext", stats)
      stats <- gsub("env_sd", "ext_sd", stats)
      return(stats)
    })
    
    cat("Loaded", length(combos), "summary stat combinations\n")
    return(combos)
  }
  
  # Legacy mode intentionally removed: correlation-based combo generation
  if (summary_stats_arg == "all") {
    stop("'all' is no longer supported. Pass an explicit combo file path instead.")
  }
  
  # Otherwise, parse as comma-separated single combination
  stats <- strsplit(summary_stats_arg, ",")[[1]]
  return(list(stats))
}

# Fixed parameters matching sample structure - load from environment or use defaults
# Defaults match the original script: 2 populations, 6 individuals, 3 replicates
num_pop <- as.integer(Sys.getenv("SIM_NUM_POP", "2"))
num_ind <- as.integer(Sys.getenv("SIM_NUM_IND", "6"))
num_rep <- as.integer(Sys.getenv("SIM_NUM_REP", "3"))

cat("Sample structure configuration:\n")
cat("  Populations:", num_pop, "\n")
cat("  Individuals per pop:", num_ind, "\n")
cat("  Replicates per ind:", num_rep, "\n")

mean_trait <- 100
var_additive <- 1

# Pre-compute indices for variance calculation (MEMORY EFFICIENT)
# These are computed once at script load, not per simulation
.pop_idx <- rep(1:num_pop, each = num_ind * num_rep)
.line_idx <- rep(1:(num_pop * num_ind), each = num_rep)
.n_total <- num_pop * num_ind * num_rep
.n_lines <- num_pop * num_ind

# Pre-compute degrees of freedom
.DF_among <- num_pop - 1
.DF_within <- num_pop * (num_ind - 1)
.DF_residual <- num_pop * num_ind * (num_rep - 1)

# Variance floor on among- and within-population ANOVA components.
# ridge_floor (default): sqrt(max(raw, 0)^2 + lambda^2), lambda = max(alpha * ANOVA noise,
#   a tiny fraction of MS_total). Components stay strictly positive. FLOOR_ALPHA default 0.1.
# baseline: clip non-positive components to that tiny MS_total fraction.
.FLOOR_POLICY <- Sys.getenv("FLOOR_POLICY", "ridge_floor")
.FLOOR_ALPHA <- as.numeric(Sys.getenv("FLOOR_ALPHA", "0.1"))
cat("  Floor policy:", .FLOOR_POLICY, " alpha=", .FLOOR_ALPHA, "\n")

.apply_floor_components <- function(var_among_raw, var_within_raw, MS_total,
                                    noise_among, noise_within,
                                    policy = .FLOOR_POLICY, alpha = .FLOOR_ALPHA) {
  if (identical(policy, "F3")) policy <- "ridge_floor"
  tf <- if (length(MS_total) == 1L) {
    if (MS_total > 0) 1e-8 * MS_total else 1e-8
  } else {
    ifelse(MS_total > 0, 1e-8 * MS_total, 1e-8)
  }
  both_neg <- (var_among_raw <= 0) & (var_within_raw <= 0)
  either_neg <- (var_among_raw <= 0) | (var_within_raw <= 0)
  if (policy == "F1") {
    lam_a <- pmax(alpha * noise_among, tf)
    lam_w <- pmax(alpha * noise_within, tf)
    var_among <- ifelse(var_among_raw > 0, var_among_raw, lam_a)
    var_within <- ifelse(var_within_raw > 0, var_within_raw, lam_w)
    flag_neg <- both_neg
  } else if (policy == "ridge_floor") {
    lam_a <- pmax(alpha * noise_among, tf)
    lam_w <- pmax(alpha * noise_within, tf)
    var_among <- sqrt(pmax(var_among_raw, 0)^2 + lam_a^2)
    var_within <- sqrt(pmax(var_within_raw, 0)^2 + lam_w^2)
    flag_neg <- rep(FALSE, length(var_among_raw))
  } else if (policy == "F4") {
    var_among <- ifelse(var_among_raw > 0, var_among_raw, tf)
    var_within <- ifelse(var_within_raw > 0, var_within_raw, tf)
    flag_neg <- either_neg
  } else if (policy == "baseline") {
    var_among <- ifelse(var_among_raw > 0, var_among_raw, tf)
    var_within <- ifelse(var_within_raw > 0, var_within_raw, tf)
    flag_neg <- both_neg
  } else {
    stop("Unknown FLOOR_POLICY='", policy, "'. Use ridge_floor (default) or baseline.")
  }
  list(var_among = var_among, var_within = var_within, both_var_negative = flag_neg)
}

#' Calculate variance components using vectorized operations (NO data.frame)
#' This is much more memory efficient than the tapply-based version
calc_variance_fast <- function(trait_values) {
  # Calculate means using matrix operations instead of tapply
  # Pop means: average of all values in each population
  # Use rowSums on matrix reshape for arbitrary num_pop
  pop_sums <- rowSums(matrix(trait_values, nrow = num_pop, ncol = num_ind * num_rep, byrow = TRUE))
  pop_means <- pop_sums / (num_ind * num_rep)
  overall_mean <- sum(trait_values) / .n_total
  
  # Line means using vectorized reshape
  line_sums <- rowSums(matrix(trait_values, nrow = .n_lines, ncol = num_rep, byrow = TRUE))
  line_means <- line_sums / num_rep
  
  # Expand means for SS calculations
  pop_means_expanded <- rep(pop_means, each = num_ind)
  line_means_expanded <- rep(line_means, each = num_rep)
  
  # Sum of squares
  SS_among <- sum((pop_means - overall_mean)^2) * num_ind * num_rep
  SS_within <- sum((line_means - pop_means_expanded)^2) * num_rep
  SS_residual <- sum((trait_values - line_means_expanded)^2)
  
  # Total SS/MS of simulated trait values (grand mean); scales with trait magnitude.
  SS_total <- sum((trait_values - overall_mean)^2)
  MS_total <- SS_total / max(.n_total - 1L, 1L)
  
  # Mean squares and variance components
  MS_among <- SS_among / .DF_among
  MS_within <- SS_within / .DF_within
  MS_residual <- SS_residual / .DF_residual
  
  var_among_raw <- (MS_among - MS_within) / (num_ind * num_rep)
  var_within_raw <- (MS_within - MS_residual) / num_rep
  var_residual <- max(MS_residual, 0)

  noise_among <- MS_within / (num_ind * num_rep)
  noise_within <- MS_residual / num_rep
  floored <- .apply_floor_components(
    var_among_raw, var_within_raw, MS_total, noise_among, noise_within
  )
  var_among <- floored$var_among
  var_within <- floored$var_within
  both_var_negative <- floored$both_var_negative

  return(c(var_among, var_within, var_residual, both_var_negative))
}

#' Generate simulated trait data and return summary statistics
#' MEMORY EFFICIENT: Uses vectorized operations, no data.frame creation
#' SELECTIVE: Only calculates the summary statistics specified in required_stats
#' 
#' @param required_stats Character vector of stats to calculate (default: all)
generate_sim_data_summarystats <- function(
  num_pop = 2, num_ind = 6, num_rep = 3,
  mean_trait = 100, sd_between_pop = 20, sd_within_pop = 10, sd_ext = 1,
  required_stats = NULL
) {
  # Default to all stats if not specified
  if (is.null(required_stats)) {
    required_stats <- ALL_SUMMARY_STATS
  }
  
  # Generate trait values directly as a vector (no data.frame)
  if (num_pop == 2) {
    mu_between_pop <- sqrt(2 * sd_between_pop^2)
    means_pop <- c(mean_trait, mean_trait + mu_between_pop)
  } else {
    means_pop <- rnorm(num_pop, mean = mean_trait, sd = sd_between_pop)
  }
  means_ind <- rep(means_pop, each = num_ind) + rnorm(num_pop * num_ind, 0, sd_within_pop)
  trait_values <- rep(means_ind, each = num_rep) + rnorm(.n_total, 0, sd_ext)
  
  # Calculate variance using fast vectorized method (always needed)
  var_components <- calc_variance_fast(trait_values)
  among_pop_var <- var_components[1]
  within_pop_var <- var_components[2]
  ext_var <- var_components[3]
  both_var_negative <- var_components[4]
  
  # Initialize result vector with only required stats
  result <- numeric(length(required_stats))
  names(result) <- required_stats
  
  # Basic stats (always calculated as part of variance components)
  if ("among_pop_sd" %in% required_stats) result["among_pop_sd"] <- sqrt(among_pop_var)
  if ("within_pop_sd" %in% required_stats) result["within_pop_sd"] <- sqrt(within_pop_var)
  if ("ext_sd" %in% required_stats) result["ext_sd"] <- sqrt(ext_var)
  if ("both_var_negative" %in% required_stats) result["both_var_negative"] <- both_var_negative
  
  # Derived stats (only calculate if needed)
  # total_genetic_var is intentionally the QST-scaled additive genetic variance:
  #   V_G = V_A_star = (V_GB + 2 * V_GW) / 2
  total_genetic_var <- (among_pop_var + 2 * within_pop_var) / 2

  if ("QST" %in% required_stats) {
    qst_denom <- among_pop_var + 2 * within_pop_var
    result["QST"] <- if (qst_denom == 0) 0 else among_pop_var / qst_denom
  }

  if ("ratioVext" %in% required_stats) {
    result["ratioVext"] <- if (total_genetic_var == 0) 0 else ext_var / total_genetic_var
  }

  if ("F_among_pop" %in% required_stats) {
    result["F_among_pop"] <- if (total_genetic_var == 0 || within_pop_var == 0) 0 else among_pop_var / within_pop_var
  }

  if ("F_within_pop" %in% required_stats) {
    result["F_within_pop"] <- if (total_genetic_var == 0 || ext_var == 0) 0 else within_pop_var / ext_var
  }

  # total_var is intentionally the QST-scaled genetic-plus-extrinsic variance:
  total_var <- total_genetic_var + ext_var
  if ("ratioVbetweenVext" %in% required_stats) {
    result["ratioVbetweenVext"] <- if (ext_var == 0) 0 else among_pop_var / ext_var
  }
  if ("ratioVbetweenVtotal" %in% required_stats) {
    result["ratioVbetweenVtotal"] <- if (total_var == 0) 0 else among_pop_var / total_var
  }
  if ("ratioVwithinVtotal" %in% required_stats) {
    result["ratioVwithinVtotal"] <- if (total_var == 0) 0 else within_pop_var / total_var
  }
  if ("ratioVextVtotal" %in% required_stats) {
    result["ratioVextVtotal"] <- if (total_var == 0) 0 else ext_var / total_var
  }

  if ("skewness_data" %in% required_stats) {
    load_e1071_if_needed()
    result["skewness_data"] <- skewness(trait_values)
  }
  if ("kurtosis_data" %in% required_stats) {
    load_e1071_if_needed()
    result["kurtosis_data"] <- kurtosis(trait_values)
  }

  return(result)
}

#' Run batch simulations - MEMORY OPTIMIZED
#' @param required_stats Character vector of stats to calculate (passed to generate_sim_data_summarystats)
run_batch_simulations <- function(batch_size, sd_between_pop_batch, sd_within_pop_batch, sd_ext_batch, required_stats = NULL) {
  if (is.null(required_stats)) required_stats <- ALL_SUMMARY_STATS
  
  # Generate traits for all simulations
  if (num_pop == 2) {
    mu_between_pop <- sqrt(2 * sd_between_pop_batch^2)
    means_pop <- cbind(rep(mean_trait, batch_size), mean_trait + mu_between_pop)
  } else {
    noise_pop <- matrix(rnorm(batch_size * num_pop), nrow = batch_size, ncol = num_pop)
    means_pop <- mean_trait + noise_pop * sd_between_pop_batch
  }
  means_pop_exp <- means_pop[, rep(1:num_pop, each = num_ind), drop = FALSE]
  
  noise_ind <- matrix(rnorm(batch_size * .n_lines), nrow = batch_size, ncol = .n_lines)
  noise_ind <- noise_ind * sd_within_pop_batch
  means_ind <- means_pop_exp + noise_ind
  
  means_ind_exp <- means_ind[, rep(1:.n_lines, each = num_rep), drop = FALSE]
  noise_rep <- matrix(rnorm(batch_size * .n_total), nrow = batch_size, ncol = .n_total)
  noise_rep <- noise_rep * sd_ext_batch
  trait_matrix <- means_ind_exp + noise_rep
  
  # ANOVA
  overall_means <- rowSums(trait_matrix) / .n_total
  
  line_sums <- matrix(0, nrow = batch_size, ncol = .n_lines)
  for (r in seq_len(num_rep)) {
    line_sums <- line_sums + trait_matrix[, seq(r, .n_total, by = num_rep), drop = FALSE]
  }
  line_means <- line_sums / num_rep
  
  pop_sums <- matrix(0, nrow = batch_size, ncol = num_pop)
  for (p in seq_len(num_pop)) {
    cols <- (p - 1) * num_ind * num_rep + seq_len(num_ind * num_rep)
    pop_sums[, p] <- rowSums(trait_matrix[, cols, drop = FALSE])
  }
  pop_means <- pop_sums / (num_ind * num_rep)
  
  SS_among <- rowSums((pop_means - overall_means)^2) * (num_ind * num_rep)
  
  SS_within <- 0
  for (p in seq_len(num_pop)) {
    line_cols <- (p - 1) * num_ind + seq_len(num_ind)
    SS_within <- SS_within + rowSums((line_means[, line_cols, drop = FALSE] - pop_means[, p])^2)
  }
  SS_within <- SS_within * num_rep
  
  SS_residual <- 0
  for (r in seq_len(num_rep)) {
    SS_residual <- SS_residual + rowSums((trait_matrix[, seq(r, .n_total, by = num_rep), drop = FALSE] - line_means)^2)
  }
  
  SS_total <- rowSums((trait_matrix - overall_means)^2)
  MS_total <- SS_total / max(.n_total - 1L, 1L)
  
  MS_among <- SS_among / .DF_among
  MS_within <- SS_within / .DF_within
  MS_residual <- SS_residual / .DF_residual
  
  var_among_raw <- (MS_among - MS_within) / (num_ind * num_rep)
  var_within_raw <- (MS_within - MS_residual) / num_rep
  var_residual <- pmax(MS_residual, 0)

  noise_among <- MS_within / (num_ind * num_rep)
  noise_within <- MS_residual / num_rep
  floored <- .apply_floor_components(
    var_among_raw, var_within_raw, MS_total, noise_among, noise_within
  )
  among_pop_var <- floored$var_among
  within_pop_var <- floored$var_within
  ext_var <- var_residual
  both_var_negative <- floored$both_var_negative

  # Compute Required Stats
  result <- matrix(0, nrow = batch_size, ncol = length(required_stats))
  colnames(result) <- required_stats
  
  if ("among_pop_sd" %in% required_stats) result[, "among_pop_sd"] <- sqrt(among_pop_var)
  if ("within_pop_sd" %in% required_stats) result[, "within_pop_sd"] <- sqrt(within_pop_var)
  if ("ext_sd" %in% required_stats) result[, "ext_sd"] <- sqrt(ext_var)
  if ("both_var_negative" %in% required_stats) result[, "both_var_negative"] <- both_var_negative
  
  total_genetic_var <- (among_pop_var + 2 * within_pop_var) / 2
  
  if ("QST" %in% required_stats) {
    qst_denom <- among_pop_var + 2 * within_pop_var
    result[, "QST"] <- ifelse(qst_denom == 0, 0, among_pop_var / qst_denom)
  }
  if ("ratioVext" %in% required_stats) {
    result[, "ratioVext"] <- ifelse(total_genetic_var == 0, 0, ext_var / total_genetic_var)
  }
  if ("F_among_pop" %in% required_stats) {
    result[, "F_among_pop"] <- ifelse(total_genetic_var == 0 | within_pop_var == 0, 0, among_pop_var / within_pop_var)
  }
  if ("F_within_pop" %in% required_stats) {
    result[, "F_within_pop"] <- ifelse(total_genetic_var == 0 | ext_var == 0, 0, within_pop_var / ext_var)
  }
  
  total_var <- total_genetic_var + ext_var
  if ("ratioVbetweenVext" %in% required_stats) {
    result[, "ratioVbetweenVext"] <- ifelse(ext_var == 0, 0, among_pop_var / ext_var)
  }
  if ("ratioVbetweenVtotal" %in% required_stats) {
    result[, "ratioVbetweenVtotal"] <- ifelse(total_var == 0, 0, among_pop_var / total_var)
  }
  if ("ratioVwithinVtotal" %in% required_stats) {
    result[, "ratioVwithinVtotal"] <- ifelse(total_var == 0, 0, within_pop_var / total_var)
  }
  if ("ratioVextVtotal" %in% required_stats) {
    result[, "ratioVextVtotal"] <- ifelse(total_var == 0, 0, ext_var / total_var)
  }
  
  if ("skewness_data" %in% required_stats) {
    load_e1071_if_needed()
    result[, "skewness_data"] <- apply(trait_matrix, 1, skewness)
  }
  if ("kurtosis_data" %in% required_stats) {
    load_e1071_if_needed()
    result[, "kurtosis_data"] <- apply(trait_matrix, 1, kurtosis)
  }
  
  return(result)
}

#' Estimate QST using ABC with MULTIPLE summary stat combinations
#' Runs simulations once and applies ABC with each combination
#' OPTIMIZED: Only calculates summary statistics that are needed
#' MEMORY EFFICIENT: Aggressive cleanup after each ABC run
#' 
#' @param obs_stats Named vector of observed summary statistics
#' @param num_sim Number of ABC simulations (default: 100000)
#' @param summary_stat_combos List of character vectors, each specifying a combination
#' @return Named list with QST estimates for each combination
estimate_qst_abc_multi <- function(obs_stats, num_sim = num_sim, summary_stat_combos, last_abc_env = NULL) {
  
  # Short-circuit: both ANOVA components non-positive -> QST not estimable (NA).
  if (is_both_neg(obs_stats)) {
    combo_names <- sapply(summary_stat_combos, function(x) paste(x, collapse = ","))
    res <- as.list(rep(NA_real_, length(combo_names)))
    names(res) <- combo_names
    return(res)
  }

  abc_seed <- abc_seed_from_obs_stats(obs_stats)
  set.seed(abc_seed)
  cat("ABC estimation seed:", abc_seed, "\n")

  # Ensure cleanup happens even on error
  
  # OPTIMIZATION: Determine which stats are actually needed
  required_stats <- get_required_stats(summary_stat_combos)
  n_stats <- length(required_stats)
  cat("Required stats (", n_stats, "/", length(ALL_SUMMARY_STATS), "):",
      paste(required_stats, collapse = ", "), "\n")
  
  batch_size <- as.integer(Sys.getenv("QST_ABC_BATCH_SIZE", "50000"))
  num_batches <- ceiling(num_sim / batch_size)
  
  # Generate prior parameters based on observed stats.
  # prior_floor prevents degenerate priors when observed genetic SDs are tiny.
  # Uses 0.1*(among + within + ext): ext_sd coupling is 0.1x, not the old 0.2x
  # that artificially inflated neutral Q_ST at high V_E/V_G.
  prior_floor <- 0.1 * (obs_stats['among_pop_sd'] + obs_stats['within_pop_sd'])
  sd_between_pop_prior <- runif(num_sim, 0, max(2 * obs_stats['among_pop_sd'], prior_floor))
  sd_within_pop_prior  <- runif(num_sim, 0, max(2 * obs_stats['within_pop_sd'], prior_floor))
  sd_ext_prior         <- runif(num_sim, 0, 2 * obs_stats['ext_sd'])
  
  # Pre-allocate result matrix - ONLY for required stats (memory optimization)
  sim_stats_matrix <- matrix(NA, nrow = num_sim, ncol = n_stats)
  prior_params <- cbind(sd_between_pop_prior, sd_within_pop_prior, sd_ext_prior)
  colnames(prior_params) <- c("sd_between_pop", "sd_within_pop", "sd_ext")
  
  rm(sd_between_pop_prior, sd_within_pop_prior, sd_ext_prior)
  
  # Process in batches
  cat("Simulating", num_sim, "for", length(summary_stat_combos), "combinations...\n")
  for (b in 1:num_batches) {
    start_idx <- (b - 1) * batch_size + 1
    end_idx <- min(b * batch_size, num_sim)
    
    batch_results <- run_batch_simulations(
      end_idx - start_idx + 1,
      prior_params[start_idx:end_idx, 1],
      prior_params[start_idx:end_idx, 2],
      prior_params[start_idx:end_idx, 3],
      required_stats = required_stats  # Pass required stats
    )
    
    sim_stats_matrix[start_idx:end_idx, ] <- batch_results
    rm(batch_results)
  }
  
  colnames(sim_stats_matrix) <- required_stats
  
  # Clean up invalid rows
  valid_rows <- complete.cases(sim_stats_matrix) & 
                (rowSums(is.infinite(sim_stats_matrix)) == 0)
  
  sim_stats_valid <- sim_stats_matrix[valid_rows, , drop = FALSE]
  prior_params_valid <- prior_params[valid_rows, , drop = FALSE]
  n_valid <- nrow(sim_stats_valid)
  
  rm(sim_stats_matrix, prior_params, valid_rows)
  
  cat("Valid simulations:", n_valid, "/", num_sim, "\n")
  
  # Prepare obs_stats with only required stats for ABC
  obs_stats_subset <- obs_stats[required_stats]
  
  # Run ABC for EACH combination on the SAME simulated data
  # MEMORY OPTIMIZATION: Aggressive cleanup after each combo (package abc path only)
  tol <- max(0.001, get_abc_tol_numerator() / n_valid)
  results <- list()
  meth_abc <- get_abc_method()

  if (meth_abc %in% c("loclinear", "rejection")) {
    cat("Vectorized fast ABC (method = ", meth_abc, ") for ", length(summary_stat_combos),
        " combinations\n", sep = "")
    vec <- fast_abc_estimate_multi(
      obs_stats_subset,
      prior_params_valid,
      sim_stats_valid,
      summary_stat_combos,
      tol = tol,
      estimator = meth_abc,
      hcorr = TRUE
    )
    results <- as.list(fast_abc_multi_qst(vec))
    if (!is.null(last_abc_env) && is.environment(last_abc_env)) {
      assign("res_abc", list(method = meth_abc, fast_vectorized = TRUE, estimates = vec),
             envir = last_abc_env)
      assign("last_error", NULL, envir = last_abc_env)
    }
  } else {
    for (i in seq_along(summary_stat_combos)) {
      combo <- summary_stat_combos[[i]]
      combo_name <- paste(combo, collapse = ",")

      qst_estimate <- tryCatch({
        combo_sumstat <- sim_stats_valid[, combo, drop = FALSE]
        combo_target <- obs_stats_subset[combo]

        res_abc <- run_abc_qst(
          target = combo_target,
          param = prior_params_valid,
          sumstat = combo_sumstat,
          tol = tol,
          transf = rep("none", ncol(prior_params_valid))
        )
        if (!is.null(last_abc_env) && is.environment(last_abc_env)) {
          assign("res_abc", res_abc, envir = last_abc_env)
          assign("last_error", NULL, envir = last_abc_env)
        }
        result <- qst_mean_from_abc(res_abc)
        rm(res_abc, combo_sumstat, combo_target)
        result

      }, error = function(e) {
        if (!is.null(last_abc_env) && is.environment(last_abc_env)) {
          assign("res_abc", NULL, envir = last_abc_env)
          assign("last_error", conditionMessage(e), envir = last_abc_env)
        }
        warning(paste("ABC failed for", combo_name, ":", e$message))
        NA
      })

      results[[combo_name]] <- qst_estimate

      # Progress indicator
      if (i %% 10 == 0 || i == length(summary_stat_combos)) {
        cat("  ABC combo", i, "/", length(summary_stat_combos), "\n")
      }
    }
  }

  rm(sim_stats_valid, prior_params_valid, obs_stats_subset)
  
  return(results)
}

#' Deterministic RNG seed for paired ABC runs (same obs_stats -> same seed).
abc_seed_from_obs_stats <- function(obs_stats) {
  env_seed <- Sys.getenv("QST_ABC_SEED", unset = "")
  if (nzchar(env_seed)) {
    v <- suppressWarnings(as.integer(env_seed))
    if (length(v) == 1L && !is.na(v)) return(v)
  }
  keys <- c("among_pop_sd", "within_pop_sd", "ext_sd")
  vals <- suppressWarnings(as.numeric(obs_stats[keys]))
  if (length(vals) != 3L || !all(is.finite(vals))) return(42L)
  raw <- paste(sprintf("%.8e", vals), collapse = "|")
  as.integer((sum(utf8ToInt(raw)) %% (.Machine$integer.max - 1L)) + 1L)
}

#' Estimate QST using ABC - disk-efficient version
#' Runs 100,000 simulations in memory, returns only QST estimate
#' 
#' @param obs_stats Named vector of observed summary statistics
#' @param num_sim Number of ABC simulations (default: 100000)
#' @param summary_stat_names Vector of summary statistic names to use for ABC
#' @return Estimated QST value
estimate_qst_abc <- function(obs_stats, num_sim = num_sim, summary_stat_names = c("QST", "ratioVext"), last_abc_env = NULL) {
  
  # Short-circuit: both ANOVA components non-positive -> QST not estimable (NA).
  if (is_both_neg(obs_stats)) {
    return(abc_na_result())
  }

  abc_seed <- abc_seed_from_obs_stats(obs_stats)
  set.seed(abc_seed)
  cat("ABC estimation seed:", abc_seed, "\n")

  # Ensure cleanup happens even on error - use reset=TRUE to return memory to OS
  
  batch_size <- as.integer(Sys.getenv("QST_ABC_BATCH_SIZE", "50000"))
  num_batches <- ceiling(num_sim / batch_size)
  
  # Generate prior parameters based on observed stats.
  # prior_floor prevents degenerate priors when observed genetic SDs are tiny.
  # Uses 0.1*(among + within + ext): ext_sd coupling is 0.1x, not the old 0.2x
  # that artificially inflated neutral Q_ST at high V_E/V_G.
  prior_floor <- 0.1 * (obs_stats['among_pop_sd'] + obs_stats['within_pop_sd'])
  sd_between_pop_prior <- runif(num_sim, 0, max(2 * obs_stats['among_pop_sd'], prior_floor))
  sd_within_pop_prior  <- runif(num_sim, 0, max(2 * obs_stats['within_pop_sd'], prior_floor))
  sd_ext_prior         <- runif(num_sim, 0, 2 * obs_stats['ext_sd'])
  
  # Pre-allocate result matrices
  # Dimensions determined by ALL_SUMMARY_STATS (defined at top of script)
  sim_stats_matrix <- matrix(NA, nrow = num_sim, ncol = length(ALL_SUMMARY_STATS))
  prior_params <- cbind(sd_between_pop_prior, sd_within_pop_prior, sd_ext_prior)
  colnames(prior_params) <- c("sd_between_pop", "sd_within_pop", "sd_ext")
  
  # Free the individual prior vectors early - they're now in prior_params
  rm(sd_between_pop_prior, sd_within_pop_prior, sd_ext_prior)
  
  # Process in batches
  cat("Simulating", num_sim, "in", num_batches, "batches...\n")
  for (b in 1:num_batches) {
    start_idx <- (b - 1) * batch_size + 1
    end_idx <- min(b * batch_size, num_sim)
    
    batch_results <- run_batch_simulations(
      end_idx - start_idx + 1,
      prior_params[start_idx:end_idx, 1],
      prior_params[start_idx:end_idx, 2],
      prior_params[start_idx:end_idx, 3]
    )
    
    sim_stats_matrix[start_idx:end_idx, ] <- batch_results
    rm(batch_results)
    # Use reset=TRUE to help return memory to OS after forked processes
    
    if (b %% 5 == 0 || b == num_batches) {
      cat("  Batch", b, "/", num_batches, "complete\n")
    }
  }
  
  colnames(sim_stats_matrix) <- ALL_SUMMARY_STATS
  
  # Clean up invalid rows - use vectorized check (faster and less memory)
  valid_rows <- complete.cases(sim_stats_matrix) & 
                (rowSums(is.infinite(sim_stats_matrix)) == 0)
  
  sim_stats_valid <- sim_stats_matrix[valid_rows, , drop = FALSE]
  prior_params_valid <- prior_params[valid_rows, , drop = FALSE]
  n_valid <- nrow(sim_stats_valid)
  
  # Free original matrices immediately after subsetting
  rm(sim_stats_matrix, prior_params, valid_rows)
  
  cat("Valid simulations:", n_valid, "/", num_sim, "\n")
  
  # Run ABC
  tol <- max(0.001, get_abc_tol_numerator() / n_valid)
  
  qst_estimate <- tryCatch({
    res_abc <- run_abc_qst(
      target = obs_stats[summary_stat_names],
      param = prior_params_valid,
      sumstat = sim_stats_valid[, summary_stat_names],
      tol = tol,
      transf = rep("none", ncol(prior_params_valid))
    )
    if (!is.null(last_abc_env) && is.environment(last_abc_env)) {
      assign("res_abc", res_abc, envir = last_abc_env)
      assign("last_error", NULL, envir = last_abc_env)
    }
    sds <- abc_sd_means_from_abc(res_abc)
    result <- abc_result_to_estimate(res_abc)
    rm(res_abc)
    result
    
  }, error = function(e) {
    if (!is.null(last_abc_env) && is.environment(last_abc_env)) {
      assign("res_abc", NULL, envir = last_abc_env)
      assign("last_error", conditionMessage(e), envir = last_abc_env)
    }
    warning(paste("ABC failed:", e$message))
    abc_na_result()
  })
  
  # Final cleanup of remaining objects
  rm(sim_stats_valid, prior_params_valid)
  # Use reset=TRUE to return memory to OS (important for forked processes in cgroups)
  
  return(qst_estimate)
}

is_valid_generated_obs <- function(obs) {
  if (is.null(obs) || !length(obs)) return(FALSE)
  if (any(is.nan(obs)) || any(is.infinite(obs))) return(FALSE)
  if ("among_pop_sd" %in% names(obs)) {
    ap <- obs["among_pop_sd"]
    if (!is.finite(ap) || ap == 0) return(FALSE)
  }
  TRUE
}

#' Generate observed summary statistics for neutral simulation
#' Uses FST to derive sd_between_pop and sd_within_pop
#' OPTIMIZED: Only calculates required stats
#' 
#' @param fst_value FST value to simulate
#' @param ratioVext V_E / V_G ratio observed from real traits
#' @param seed Random seed for reproducibility
#' @param required_stats Character vector of stats to calculate (default: all)
generate_neutral_obs_stats <- function(fst_value, ratioVext, seed = NULL, required_stats = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Derive variance parameters from FST
  sd_within_pop <- sqrt(abs((1 - fst_value) * var_additive))
  sd_between_pop <- sqrt(abs(2 * fst_value * var_additive))
  
  # Derive ext_sd properly mapped via V_E / V_G ratio to match simulated additive variance limits
  ext_sd <- sqrt(ratioVext * var_additive)
  
  # Generate "observed" data for this neutral scenario
  obs_valid <- FALSE
  attempts <- 0
  max_attempts <- 100
  
  while (!obs_valid && attempts < max_attempts) {
    obs <- generate_sim_data_summarystats(
      num_pop = num_pop, num_ind = num_ind, num_rep = num_rep,
      mean_trait = mean_trait,
      sd_between_pop = sd_between_pop,
      sd_within_pop = sd_within_pop,
      sd_ext = ext_sd,
      required_stats = required_stats
    )
    
    if (is_valid_generated_obs(obs)) {
      obs_valid <- TRUE
    }
    attempts <- attempts + 1
  }
  
  if (!obs_valid) {
    warning("Could not generate valid neutral observation")
    return(NULL)
  }
  
  return(obs)
}

#' Generate observed summary statistics for adaptive trait simulation (Design Module)
#' Uses pre-defined QST and V_E/V_G ratio to derive variance parameters
#' OPTIMIZED: Only calculates required stats
#' 
#' @param qst_value Pre-defined QST value (0-1)
#' @param ve_ratio V_E / V_G ratio (extrinsic variance / genetic variance)
#' @param seed Random seed for reproducibility
#' @param required_stats Character vector of stats to calculate (default: all)
#' @return Named vector of summary statistics, or NULL if generation fails
generate_adaptive_obs_stats <- function(qst_value, ve_ratio, seed = NULL, required_stats = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Derive variance parameters from QST
  # QST = V_GB / (V_GB + 2*V_GW), so if we fix total genetic variance = var_additive
  # V_GB = QST * (V_GB + 2*V_GW) => V_GB = QST * V_total_genetic / (1 - QST + 2*QST) 
  # Simpler: use same formula as neutral (FST = QST under neutrality)
  sd_within_pop <- sqrt(abs((1 - qst_value) * var_additive))
  sd_between_pop <- sqrt(abs(2 * qst_value * var_additive))
  
  # Derive ext_sd from V_E ratio
  # ve_ratio = V_E / V_G, where V_G = V_GB + V_GW = var_additive
  ext_sd <- sqrt(ve_ratio * var_additive)
  
  # Generate "observed" data for this adaptive trait
  obs_valid <- FALSE
  attempts <- 0
  max_attempts <- 100
  
  while (!obs_valid && attempts < max_attempts) {
    obs <- generate_sim_data_summarystats(
      num_pop = num_pop, num_ind = num_ind, num_rep = num_rep,
      mean_trait = mean_trait,
      sd_between_pop = sd_between_pop,
      sd_within_pop = sd_within_pop,
      sd_ext = ext_sd,
      required_stats = required_stats
    )
    
    if (is_valid_generated_obs(obs)) {
      obs_valid <- TRUE
    }
    attempts <- attempts + 1
  }
  
  if (!obs_valid) {
    warning(paste("Could not generate valid adaptive observation for QST =", qst_value))
    return(NULL)
  }
  
  return(obs)
}

#' Direct ANOVA QST estimate from observed summary statistics (Design Module perf-eval)
estimate_qst_anova <- function(obs_stats) {
  if (is_both_neg(obs_stats)) return(NA_real_)
  if (!"QST" %in% names(obs_stats)) return(NA_real_)
  qst <- suppressWarnings(as.numeric(obs_stats["QST"]))
  if (length(qst) != 1L || is.na(qst)) return(NA_real_)
  qst
}

#' 2x6 without-replication variance components (soft-clip floor, no residual df)
calc_variance_norep_fast <- function(trait_values) {
  pop_sums <- rowSums(matrix(trait_values, nrow = num_pop, ncol = num_ind, byrow = TRUE))
  pop_means <- pop_sums / num_ind
  overall_mean <- mean(trait_values)

  SS_among <- sum((pop_means - overall_mean)^2) * num_ind
  pop_means_expanded <- rep(pop_means, each = num_ind)
  SS_within <- sum((trait_values - pop_means_expanded)^2)

  MS_among <- SS_among / (num_pop - 1)
  MS_within <- SS_within / (num_pop * (num_ind - 1))

  var_among_raw <- (MS_among - MS_within) / num_ind
  var_within_raw <- MS_within

  SS_total <- sum((trait_values - overall_mean)^2)
  MS_total <- SS_total / max(num_pop * num_ind - 1L, 1L)
  noise_among <- MS_within / num_ind
  noise_within <- MS_within
  floored <- .apply_floor_components(
    var_among_raw, var_within_raw, MS_total, noise_among, noise_within
  )
  c(floored$var_among, floored$var_within, floored$both_var_negative)
}

#' Generate simulated trait summary stats without replication (2 pop x 6 strain)
generate_sim_data_summarystats_norep <- function(
  num_pop = 2, num_ind = 6,
  mean_trait = 100, sd_between_pop = 20, sd_within_pop = 10, sd_ext = 1,
  required_stats = NULL
) {
  if (is.null(required_stats)) {
    required_stats <- c("among_pop_sd", "within_pop_sd", "QST", "both_var_negative")
  }

  if (num_pop == 2) {
    mu_between_pop <- sqrt(2 * sd_between_pop^2)
    means_pop <- c(mean_trait, mean_trait + mu_between_pop)
  } else {
    means_pop <- rnorm(num_pop, mean = mean_trait, sd = sd_between_pop)
  }
  n_total <- num_pop * num_ind
  means_ind <- rep(means_pop, each = num_ind) +
    rnorm(n_total, 0, sd_within_pop) +
    rnorm(n_total, 0, sd_ext)

  var_components <- calc_variance_norep_fast(means_ind)
  among_pop_var <- var_components[1]
  within_pop_var <- var_components[2]
  both_var_negative <- var_components[3]

  result <- numeric(length(required_stats))
  names(result) <- required_stats
  if ("among_pop_sd" %in% required_stats) result["among_pop_sd"] <- sqrt(among_pop_var)
  if ("within_pop_sd" %in% required_stats) result["within_pop_sd"] <- sqrt(within_pop_var)
  if ("both_var_negative" %in% required_stats) result["both_var_negative"] <- both_var_negative
  if ("QST" %in% required_stats) {
    qst_denom <- among_pop_var + 2 * within_pop_var
    result["QST"] <- if (qst_denom == 0) 0 else among_pop_var / qst_denom
  }
  result
}

generate_neutral_obs_stats_norep <- function(fst_value, ratioVext, seed = NULL, required_stats = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sd_within_pop <- sqrt(abs((1 - fst_value) * var_additive))
  sd_between_pop <- sqrt(abs(2 * fst_value * var_additive))
  ext_sd <- sqrt(ratioVext * var_additive)

  obs_valid <- FALSE
  attempts <- 0
  max_attempts <- 100
  while (!obs_valid && attempts < max_attempts) {
    obs <- generate_sim_data_summarystats_norep(
      num_pop = num_pop, num_ind = num_ind,
      mean_trait = mean_trait,
      sd_between_pop = sd_between_pop,
      sd_within_pop = sd_within_pop,
      sd_ext = ext_sd,
      required_stats = required_stats
    )
    if (is_valid_generated_obs(obs)) {
      obs_valid <- TRUE
    }
    attempts <- attempts + 1
  }
  if (!obs_valid) {
    warning("Could not generate valid neutral observation (norep)")
    return(NULL)
  }
  obs
}

generate_adaptive_obs_stats_norep <- function(qst_value, ve_ratio, seed = NULL, required_stats = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sd_within_pop <- sqrt(abs((1 - qst_value) * var_additive))
  sd_between_pop <- sqrt(abs(2 * qst_value * var_additive))
  ext_sd <- sqrt(ve_ratio * var_additive)

  obs_valid <- FALSE
  attempts <- 0
  max_attempts <- 100
  while (!obs_valid && attempts < max_attempts) {
    obs <- generate_sim_data_summarystats_norep(
      num_pop = num_pop, num_ind = num_ind,
      mean_trait = mean_trait,
      sd_between_pop = sd_between_pop,
      sd_within_pop = sd_within_pop,
      sd_ext = ext_sd,
      required_stats = required_stats
    )
    if (is_valid_generated_obs(obs)) {
      obs_valid <- TRUE
    }
    attempts <- attempts + 1
  }
  if (!obs_valid) {
    warning(paste("Could not generate valid adaptive observation (norep) for QST =", qst_value))
    return(NULL)
  }
  obs
}

#' Process a batch of repeats for adaptive QST estimation (Design Module)
#' Simulates traits with pre-defined QST and V_E ratio, estimates QST via ABC
#' 
#' @param n_repeats Number of repeats to process
#' @param qst_value Pre-defined QST value
#' @param ve_ratio V_E / V_G ratio
#' @param num_sim Number of ABC simulations per repeat
#' @param summary_stat_names Summary statistics to use for ABC
#' @param start_id Starting ID for reproducible seeds
#' @return Data frame with repeat_id, true_qst, ve_ratio, and estimated_qst columns
process_evaluate_batch <- function(n_repeats, qst_value, ve_ratio, num_sim, 
                                   summary_stat_names = c("QST", "ratioVext"),
                                   start_id = 1) {
  results <- data.frame(
    repeat_id = seq(start_id, start_id + n_repeats - 1),
    true_qst = rep(qst_value, n_repeats),
    ve_ratio = rep(ve_ratio, n_repeats),
    estimated_qst = rep(NA_real_, n_repeats)
  )
  
  for (i in 1:n_repeats) {
    repeat_id <- results$repeat_id[i]
    
    # Set seed for reproducibility (unique per repeat/qst/ratio combination)
    seed <- (round(qst_value * 100) * 10000 + round(ve_ratio * 1000) * 100 + repeat_id) %% .Machine$integer.max
    
    # Generate adaptive observed stats
    obs_stats <- generate_adaptive_obs_stats(qst_value, ve_ratio, seed)
    
    if (!is.null(obs_stats)) {
      # Run ABC estimation
      qst_estimate <- estimate_qst_abc(obs_stats, num_sim, summary_stat_names)
      results$estimated_qst[i] <- qst_estimate$qst
    }
    
    # Force garbage collection after each repeat to prevent memory accumulation
    
    # Progress update every 10 repeats
    if (i %% 10 == 0 || i == n_repeats) {
      cat("  Processed", i, "/", n_repeats, "repeats\n")
    }
  }
  
  return(results)
}

#' Process a batch of repeats with MULTIPLE summary stat combinations (Design Module)
#' Runs ABC with all combinations on the SAME simulated data for efficiency
#' OPTIMIZED: Only calculates required summary statistics
#' 
#' @param n_repeats Number of repeats to process
#' @param qst_value Pre-defined QST value
#' @param ve_ratio V_E / V_G ratio
#' @param num_sim Number of ABC simulations per repeat
#' @param summary_stat_combos List of character vectors, each specifying a combination
#' @param start_id Starting ID for reproducible seeds
#' @return Data frame with repeat_id, true_qst, ve_ratio, combo, and estimated_qst columns
process_evaluate_batch_multi <- function(n_repeats, qst_value, ve_ratio, num_sim, 
                                         summary_stat_combos,
                                         start_id = 1) {
  n_combos <- length(summary_stat_combos)
  combo_names <- sapply(summary_stat_combos, function(x) paste(x, collapse = ","))
  
  # OPTIMIZATION: Determine required stats upfront
  required_stats <- get_required_stats(summary_stat_combos)
  cat("Required stats for batch:", paste(required_stats, collapse = ", "), "\n")
  
  # Create result data frame with all combinations
  results <- data.frame(
    repeat_id = rep(seq(start_id, start_id + n_repeats - 1), each = n_combos),
    true_qst = rep(qst_value, n_repeats * n_combos),
    ve_ratio = rep(ve_ratio, n_repeats * n_combos),
    combo = rep(combo_names, n_repeats),
    estimated_qst = rep(NA_real_, n_repeats * n_combos),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n_repeats) {
    repeat_id <- start_id + i - 1
    
    # Set seed for reproducibility
    seed <- (round(qst_value * 100) * 10000 + round(ve_ratio * 1000) * 100 + repeat_id) %% .Machine$integer.max
    
    # Generate adaptive observed stats (only required ones)
    obs_stats <- generate_adaptive_obs_stats(qst_value, ve_ratio, seed, required_stats)
    
    if (!is.null(obs_stats)) {
      # Run ABC with ALL combinations on the SAME simulated data
      qst_estimates <- estimate_qst_abc_multi(obs_stats, num_sim, summary_stat_combos)
      
      # Store results for each combination
      for (j in seq_along(combo_names)) {
        idx <- (i - 1) * n_combos + j
        results$estimated_qst[idx] <- qst_estimates[[combo_names[j]]]
      }
    }
    
    if (i %% 10 == 0 || i == n_repeats) {
      cat("  Processed", i, "/", n_repeats, "repeats\n")
    }
  }
  
  return(results)
}

#' Process a batch of FST values with MULTIPLE summary stat combinations
#' Runs ABC with all combinations on the SAME simulated data for efficiency
#' OPTIMIZED: Only calculates required summary statistics
#' 
#' @param fst_values Vector of FST values to process
#' @param ratioVext V_E / V_G ratio observed from real traits
#' @param num_sim Number of ABC simulations per FST
#' @param summary_stat_combos List of character vectors, each specifying a combination
process_fst_batch_multi <- function(fst_values, ratioVext, num_sim, summary_stat_combos, start_idx = 1) {
  fast_process_fst_batch_multi(fst_values, ratioVext, num_sim, summary_stat_combos, start_idx)
}

#' Process a batch of FST values for neutral QST estimation
#' Processes FST values sequentially (ABC already uses parallelism internally)
#' 
#' @param fst_values Vector of FST values to process
#' @param ratioVext V_E / V_G ratio observed from real traits
#' @param num_sim Number of ABC simulations per FST
#' @param summary_stat_names Summary statistics to use for ABC
#' @return Data frame with fst and qst columns
process_fst_batch <- function(fst_values, ratioVext, num_sim, summary_stat_names = c("QST", "ratioVext"), start_idx = 1) {
  fast_process_fst_batch(fst_values, ratioVext, num_sim, summary_stat_names, start_idx)
}

# Command-line interface (skip when sourcing for diagnostics: Sys.setenv(QST_ABC_SOURCE_ONLY=1))
if (!interactive() && Sys.getenv("QST_ABC_SOURCE_ONLY", "") != "1") {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 4) {
    cat("Usage: Rscript qst_abc_sim.R <mode> <input> <ratioVext_or_ve_ratio> <output_file> [num_sim] [summary_stats]\n")
    cat("\nModes:\n")
    cat("  Detection Module:\n")
    cat("    trait        - Estimate QST for observed trait data\n")
    cat("    neutral      - Estimate QST for single neutral FST value\n")
    cat("    batch_neutral - Estimate QST for batch of neutral FST values\n")
    cat("  Design Module (Performance Evaluation):\n")
    cat("    evaluate     - Estimate QST for single adaptive QST/VE_ratio\n")
    cat("    batch_evaluate - Estimate QST for batch of repeats with same QST/VE_ratio\n")
    cat("\nArguments:\n")
    cat("  input: For trait mode: path to observed stats RData\n")
    cat("         For neutral mode: FST value\n")
    cat("         For batch_neutral mode: path to FST file\n")
    cat("         For evaluate mode: QST value (e.g., 0.75)\n")
    cat("         For batch_evaluate mode: 'start_id_n_repeats_qst_value'\n")
    cat("  ratioVext_or_ve_ratio: For trait mode: ignored (read from input file)\n")
    cat("                      For neutral/batch_neutral mode: ratioVext extracted from true trait observations\n")
    cat("                      For evaluate/batch_evaluate mode: V_E/V_G ratio\n")
    cat("  output_file: Path to save result\n")
    cat("  num_sim: Number of ABC simulations (default: 100000)\n")
    cat("  summary_stats: Can be one of:\n")
    cat("                 - Comma-separated stats: 'QST,ratioVbetweenVtotal' (default)\n")
    cat("                 - Path to file with tab-separated combinations\n")
    quit(status = 1)
  }
  
  mode <- args[1]
  input <- args[2]
  ratioVext_arg <- args[3]
  output_file <- args[4]
  num_sim <- if (length(args) >= 5) as.numeric(args[5]) else 100000
  summary_stats_arg <- if (length(args) >= 6) args[6] else "QST,ratioVbetweenVtotal"
  
  # Parse summary stats argument - comma-string or file path
  summary_stat_combos <- parse_summary_stats(summary_stats_arg)
  is_multi_combo <- !is.null(summary_stat_combos) && length(summary_stat_combos) > 1
  
  # For single combination, extract the vector
  if (!is.null(summary_stat_combos) && length(summary_stat_combos) == 1) {
    summary_stats <- summary_stat_combos[[1]]
  } else {
    summary_stats <- summary_stat_combos[[1]]  # Use first for display
  }
  
  cat("Mode:", mode, "\n")
  cat("Input:", input, "\n")
  cat("Num simulations:", num_sim, "\n")
  if (is_multi_combo) {
    cat("Summary stats:", length(summary_stat_combos), "combinations\n")
  } else {
    cat("Summary stats:", paste(summary_stats, collapse = ", "), "\n")
  }
  cat("Using", num_cores, "cores\n")
  
  start_time <- Sys.time()
  
  if (mode == "trait") {
    # Load observed summary statistics
    load(input)  # Expects 'obs_stats' variable
    cat("Ext SD from obs_stats:", obs_stats['ext_sd'], "\n")
    
    both_neg <- is_both_neg(obs_stats)
    if (both_neg) {
      cat("NOTE: both ANOVA variance components non-positive -> trait QST set to NA\n")
      abc_result <- abc_na_result()
    } else {
      abc_result <- estimate_qst_abc(obs_stats, num_sim, summary_stats)
    }
    
    result <- list(
      mode = "trait",
      qst = abc_result$qst,
      qst_lo = abc_result$qst_lo,
      qst_hi = abc_result$qst_hi,
      abc_among_pop_sd = abc_result$abc_among_pop_sd,
      abc_within_pop_sd = abc_result$abc_within_pop_sd,
      abc_ext_sd = abc_result$abc_ext_sd,
      abc_ratioVext = abc_result$abc_ratioVext,
      abc_ratioVext_lo = abc_result$abc_ratioVext_lo,
      abc_ratioVext_hi = abc_result$abc_ratioVext_hi,
      anova_ratioVext = as.numeric(obs_stats["ratioVext"]),
      ratioVext_source = "abc",
      both_var_negative = as.integer(both_neg)
    )
    
  } else if (mode == "neutral") {
    fst_value <- as.numeric(input)
    
    rv_info <- resolve_ratioVext_arg(ratioVext_arg)
    ratioVext <- rv_info$ratioVext
    
    if (!isTRUE(rv_info$valid)) {
      result <- list(mode = "neutral", fst = fst_value, qst = NA)
    } else {
      # Set seed based on FST value for reproducibility
      seed <- round(fst_value * 1e6) %% .Machine$integer.max
      
      # Generate neutral observed stats
      obs_stats <- generate_neutral_obs_stats(fst_value, ratioVext, seed)
      
      if (is.null(obs_stats)) {
        result <- list(mode = "neutral", fst = fst_value, qst = NA)
      } else {
        qst_estimate <- estimate_qst_abc(obs_stats, num_sim, summary_stats)
        result <- list(
          mode = "neutral",
          fst = fst_value,
          qst = qst_estimate$qst
        )
      }
    }
    
  } else if (mode == "batch_neutral") {
    # Batch mode: process multiple FST values in one job
    # ratioVext_arg: numeric VE/VG, obs_stats path, or trait_qst.RData (abcVEVG)
    
    rv_info <- resolve_ratioVext_arg(ratioVext_arg)
    ratioVext <- rv_info$ratioVext
    
    # Parse input - either "start:end:file" or just file path
    if (grepl(":", input)) {
      parts <- strsplit(input, ":")[[1]]
      start_idx <- as.integer(parts[1])
      end_idx <- as.integer(parts[2])
      fst_file <- parts[3]
      
      # Read FST values from file
      all_fst <- scan(fst_file, what = numeric(), quiet = TRUE)
      fst_values <- all_fst[start_idx:end_idx]
      cat("Processing FST values", start_idx, "to", end_idx, "from", fst_file, "\n")
    } else {
      # Input is a file with FST values
      fst_values <- scan(input, what = numeric(), quiet = TRUE)
      cat("Processing", length(fst_values), "FST values from", input, "\n")
      
      # Extract batch number from filename to ensure unique seeds across batches
      batch_num <- suppressWarnings(as.integer(gsub(".*fst_batch_([0-9]+)\\.txt.*", "\\1", input)))
      if (is.na(batch_num)) batch_num <- 1
      start_idx <- (batch_num - 1) * 10000 + 1
    }
    
    cat("Number of FST values:", length(fst_values), "\n")
    
    if (!isTRUE(rv_info$valid)) {
      cat("NOTE: invalid ratioVext from trait ABC -> all neutral QST set to NA\n")
      batch_results <- data.frame(
        fst = fst_values,
        qst = rep(NA_real_, length(fst_values)),
        both_var_negative = rep(1L, length(fst_values))
      )
    } else if (is_multi_combo) {
      cat("Running ABC with", length(summary_stat_combos), "summary stat combinations\n")
      batch_results <- process_fst_batch_multi(fst_values, ratioVext, num_sim, summary_stat_combos, start_idx)
    } else {
      batch_results <- process_fst_batch(fst_values, ratioVext, num_sim, summary_stats, start_idx)
    }
    
    result <- list(
      mode = "batch_neutral",
      n_fst = length(fst_values),
      is_multi_combo = is_multi_combo,
      results = batch_results
    )
    
  } else if (mode == "full_trait") {
    fst_file <- ratioVext_arg
    load(input)
    cat("Ext SD from obs_stats:", obs_stats["ext_sd"], "\n")

    both_neg <- is_both_neg(obs_stats)
    if (both_neg) {
      cat("NOTE: both ANOVA variance components non-positive -> trait QST set to NA\n")
      abc_result <- abc_na_result()
    } else {
      abc_result <- estimate_qst_abc(obs_stats, num_sim, summary_stats)
    }

    trait_result <- list(
      mode = "trait",
      qst = abc_result$qst,
      qst_lo = abc_result$qst_lo,
      qst_hi = abc_result$qst_hi,
      abc_among_pop_sd = abc_result$abc_among_pop_sd,
      abc_within_pop_sd = abc_result$abc_within_pop_sd,
      abc_ext_sd = abc_result$abc_ext_sd,
      abc_ratioVext = abc_result$abc_ratioVext,
      abc_ratioVext_lo = abc_result$abc_ratioVext_lo,
      abc_ratioVext_hi = abc_result$abc_ratioVext_hi,
      anova_ratioVext = as.numeric(obs_stats["ratioVext"]),
      ratioVext_source = "abc",
      both_var_negative = as.integer(both_neg),
      abc_na = as.integer(!both_neg && is.na(abc_result$qst))
    )

    trait_qst_path <- sub("_result\\.csv$", "_trait_qst.RData", output_file)
    result <- trait_result
    save(result, file = trait_qst_path)
    cat("Saved trait ABC:", trait_qst_path, "\n")

    fst_values <- scan(fst_file, what = numeric(), quiet = TRUE)
    cat("Processing", length(fst_values), "FST values from", fst_file, "\n")
    batch_num <- suppressWarnings(as.integer(gsub(".*fst_batch_([0-9]+)\\.txt.*", "\\1", fst_file)))
    if (is.na(batch_num)) batch_num <- 1
    start_idx <- (batch_num - 1) * 10000 + 1

    ratioVext <- trait_result$abc_ratioVext
    ratio_valid <- !is.null(ratioVext) && length(ratioVext) == 1L && is.finite(ratioVext)
    if (!ratio_valid) {
      cat("NOTE: invalid abc_ratioVext from trait ABC -> all neutral QST set to NA\n")
      batch_results <- data.frame(
        fst = fst_values,
        qst = rep(NA_real_, length(fst_values)),
        both_var_negative = rep(1L, length(fst_values))
      )
    } else if (is_multi_combo) {
      batch_results <- process_fst_batch_multi(fst_values, ratioVext, num_sim, summary_stat_combos, start_idx)
    } else {
      batch_results <- process_fst_batch(fst_values, ratioVext, num_sim, summary_stats, start_idx)
    }

    neutral_path <- "neutral_batch_1.RData"
    result <- list(
      mode = "batch_neutral",
      n_fst = length(fst_values),
      is_multi_combo = is_multi_combo,
      results = batch_results
    )
    save(result, file = neutral_path)
    cat("Saved neutral batch:", neutral_path, "\n")

    agg_script <- file.path(script_dir, "aggregate_qst.R")
    if (!file.exists(agg_script)) stop("Missing aggregate script: ", agg_script)
    old_src <- Sys.getenv("QST_ABC_SOURCE_ONLY", "")
    Sys.setenv(QST_ABC_SOURCE_ONLY = "1")
    sys.source(agg_script, envir = environment(), keep.source = FALSE)
    if (nzchar(old_src)) Sys.setenv(QST_ABC_SOURCE_ONLY = old_src) else Sys.unsetenv("QST_ABC_SOURCE_ONLY")

    threshold_pct <- as.numeric(Sys.getenv("QST_THRESHOLD_PERCENTILE", "0.95"))
    sanity <- identical(toupper(Sys.getenv("SANITY_CHECK", "FALSE")), "TRUE")
    aggregate_qst_from_files(trait_qst_path, dirname(trait_qst_path), threshold_pct, output_file, sanity)
    result <- list(mode = "full_trait", output_file = output_file)

  } else if (mode == "evaluate") {
    # Design Module: Single adaptive QST evaluation
    # Input: QST value, ratioVext_arg: V_E/V_G ratio
    
    qst_value <- as.numeric(input)
    ve_ratio <- as.numeric(ratioVext_arg)
    
    cat("Adaptive QST:", qst_value, "\n")
    cat("V_E/V_G ratio:", ve_ratio, "\n")
    
    # Set seed for reproducibility
    seed <- (round(qst_value * 100) * 1000 + round(ve_ratio * 1000)) %% .Machine$integer.max
    
    # Generate adaptive observed stats
    obs_stats <- generate_adaptive_obs_stats(qst_value, ve_ratio, seed)
    
    if (is.null(obs_stats)) {
      result <- list(mode = "evaluate", true_qst = qst_value, ve_ratio = ve_ratio, 
                     estimated_qst = NA)
    } else {
      # Run ABC estimation
      qst_estimate <- estimate_qst_abc(obs_stats, num_sim, summary_stats)
      result <- list(
        mode = "evaluate",
        true_qst = qst_value,
        ve_ratio = ve_ratio,
        estimated_qst = qst_estimate$qst
      )
    }
    
  } else if (mode == "batch_evaluate") {
    # Design Module: Batch adaptive QST evaluation
    # Input format: "n_repeats_qst_value" OR "start_id_n_repeats_qst_value" (using underscore separator)
    # ratioVext_arg: V_E/V_G ratio
    
    ve_ratio <- as.numeric(ratioVext_arg)
    cat("V_E/V_G ratio:", ve_ratio, "\n")
    
    # Parse input (use underscore as separator to avoid HTCondor colon issues)
    parts <- strsplit(input, "_")[[1]]
    if (length(parts) == 2) {
      n_repeats <- as.integer(parts[1])
      qst_value <- as.numeric(parts[2])
      start_id <- 1
    } else if (length(parts) == 3) {
      start_id <- as.integer(parts[1])
      n_repeats <- as.integer(parts[2])
      qst_value <- as.numeric(parts[3])
    } else {
      stop("Invalid input format for batch_evaluate. Expected 'n_repeats_qst' or 'start_id_n_repeats_qst'")
    }
    
    cat("Adaptive QST:", qst_value, "\n")
    cat("Number of repeats:", n_repeats, "\n")
    cat("Starting ID:", start_id, "\n")
    
    # Process all repeats - use multi-combo if multiple combinations
    if (is_multi_combo) {
      cat("Running ABC with", length(summary_stat_combos), "summary stat combinations\n")
      batch_results <- process_evaluate_batch_multi(n_repeats, qst_value, ve_ratio, num_sim, 
                                                    summary_stat_combos, start_id)
    } else {
      batch_results <- process_evaluate_batch(n_repeats, qst_value, ve_ratio, num_sim, 
                                              summary_stats, start_id)
    }
    
    result <- list(
      mode = "batch_evaluate",
      true_qst = qst_value,
      ve_ratio = ve_ratio,
      n_repeats = n_repeats,
      start_id = start_id,
      is_multi_combo = is_multi_combo,
      results = batch_results
    )
    
  } else {
    stop("Unknown mode: ", mode)
  }
  
  if (mode != "full_trait") {
    save(result, file = output_file)
  }
  
  end_time <- Sys.time()
  elapsed_mins <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  if (mode == "full_trait") {
    cat("Fused full_trait complete; result CSV:", output_file, "\n")
  } else if (mode == "batch_neutral") {
    n_valid <- sum(!is.na(result$results$qst))
    cat("Batch complete:", n_valid, "/", result$n_fst, "valid QST estimates\n")
    cat("Average time per FST:", round(elapsed_mins * 60 / result$n_fst, 2), "seconds\n")
  } else if (mode == "batch_evaluate") {
    n_valid <- sum(!is.na(result$results$estimated_qst))
    cat("Batch complete:", n_valid, "/", result$n_repeats, "valid QST estimates\n")
    cat("Average time per repeat:", round(elapsed_mins * 60 / result$n_repeats, 2), "seconds\n")
    cat("True QST:", result$true_qst, "| Mean estimated QST:", 
        round(mean(result$results$estimated_qst, na.rm = TRUE), 4), "\n")
  } else if (mode == "evaluate") {
    cat("True QST:", result$true_qst, "| Estimated QST:", result$estimated_qst, "\n")
  } else {
    cat("QST estimate:", result$qst, "\n")
  }
  cat("Completed in", round(elapsed_mins, 2), "minutes\n")
  cat("Saved to:", output_file, "\n")
}
