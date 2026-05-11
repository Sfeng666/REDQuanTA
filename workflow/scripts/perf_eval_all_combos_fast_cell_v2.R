#!/usr/bin/env Rscript
# Fast all-combo perf-eval worker (v2: vectorized hot paths).
#
# Differences vs v1 (perf_eval_all_combos_fast_cell.R):
#   - Precompute qst_pool once per dataset (rejection mode no longer divides per combo).
#   - Inner combo loop works on squared distances (no sqrt per combo).
#   - loclinear_adjust uses colMeans + sweep(..., check.margin=FALSE) and avoids
#     materializing full N x K pred / pred_sd matrices (they are constant per column).
#   - Distance chunk uses tcrossprod(diff2, mask) instead of diff2 %*% t(mask).
#   - Per-dataset MAD scaling uses vapply over columns instead of apply.
#   - Output is a named numeric vector rather than a list.
# Same FAST_ABC_ESTIMATOR contract (loclinear default, rejection optional) and
# same RData schema as v1.
#
# Optional tool: not wired into the default REDQuanTA Snakemake workflow; place
# qst_abc_sim.R on QST_CODE_DIR or alongside this script when running standalone.

get_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", ca, value = TRUE)
  if (length(f)) {
    return(dirname(normalizePath(sub("^--file=", "", f[[1]]), winslash = "/")))
  }
  getwd()
}

source_core <- function() {
  code_dir <- Sys.getenv("QST_CODE_DIR", unset = "")
  if (!nzchar(code_dir)) code_dir <- get_script_dir()
  core <- file.path(code_dir, "qst_abc_sim.R")
  if (!file.exists(core)) stop("qst_abc_sim.R not found: ", core)

  old <- Sys.getenv("QST_ABC_SOURCE_ONLY", unset = NA)
  Sys.setenv(QST_ABC_SOURCE_ONLY = "1")
  source(core, local = FALSE)
  if (is.na(old)) {
    Sys.unsetenv("QST_ABC_SOURCE_ONLY")
  } else {
    Sys.setenv(QST_ABC_SOURCE_ONLY = old)
  }
}

source_core()

fast_estimator_mode <- function(estimator = NULL) {
  if (!is.null(estimator)) return(estimator)
  mode <- Sys.getenv("FAST_ABC_ESTIMATOR", unset = "loclinear")
  if (!mode %in% c("loclinear", "rejection")) {
    stop("FAST_ABC_ESTIMATOR must be loclinear or rejection, got: ", mode)
  }
  mode
}

mad_scale <- function(x) {
  s <- stats::mad(x)
  if (!is.finite(s) || s == 0) 1 else s
}

mad_scales <- function(mat) {
  vapply(seq_len(ncol(mat)), function(j) mad_scale(mat[, j]), numeric(1))
}

# Returns the indices of the smallest `nacc = ceiling(N * tol)` values of d2.
# Operates on squared distances (ordering is preserved under sqrt).
accepted_indices_sq <- function(d2, tol) {
  n <- length(d2)
  nacc <- max(1L, min(n, ceiling(n * tol)))
  ds <- sort(d2, partial = nacc)[nacc]
  wt <- d2 <= ds
  aux <- cumsum(wt)
  which(wt & aux <= nacc)
}

epanechnikov_weights_sq <- function(d2_pick) {
  d2max <- max(d2_pick)
  if (!is.finite(d2max) || d2max <= 0) return(rep(1, length(d2_pick)))
  pmax(0, 1 - d2_pick / d2max)
}

qst_from_params_vec <- function(sb_pick, sw_pick) {
  sb2 <- sb_pick * sb_pick
  mean(sb2 / (sb2 + 2 * sw_pick * sw_pick), na.rm = TRUE)
}

qst_from_param_matrix <- function(param) {
  if (is.null(param) || nrow(param) < 1L) return(NA_real_)
  qst_from_params_vec(param[, "sd_between_pop"], param[, "sd_within_pop"])
}

# Vectorized loclinear adjustment mirroring abc::abc(method = "loclinear")
# after accepted-row selection and Epanechnikov weighting.
loclinear_adjust <- function(target_scaled, param_pick, sumstat_pick_scaled, weights,
                             hcorr = TRUE) {
  X <- cbind(1, sumstat_pick_scaled)
  target_with_intercept <- c(1, target_scaled)

  fit1 <- stats::lsfit(sumstat_pick_scaled, param_pick, wt = weights)
  coef1 <- cbind(fit1$coefficients)[fit1$qr$pivot, , drop = FALSE]
  pred_at_target <- as.numeric(crossprod(coef1, target_with_intercept))

  fitted_vals <- X %*% coef1
  residuals <- param_pick - fitted_vals
  resid_mean <- colMeans(residuals)
  residuals <- sweep(residuals, 2, resid_mean, "-", check.margin = FALSE)
  pred_at_target <- pred_at_target + resid_mean

  if (isTRUE(hcorr)) {
    corrected <- tryCatch({
      log_resid2 <- log(pmax(residuals * residuals, .Machine$double.xmin))
      fit2 <- stats::lsfit(sumstat_pick_scaled, log_resid2, wt = weights)
      coef2 <- cbind(fit2$coefficients)[fit2$qr$pivot, , drop = FALSE]
      log_pred_sd_target <- as.numeric(crossprod(coef2, target_with_intercept))
      pred_sd_vec <- sqrt(exp(log_pred_sd_target))
      log_pred_si <- X %*% coef2
      pred_si <- sqrt(exp(log_pred_si))
      tmp <- residuals / pred_si
      tmp <- sweep(tmp, 2, pred_sd_vec, "*", check.margin = FALSE)
      adj <- sweep(tmp, 2, pred_at_target, "+", check.margin = FALSE)
      if (any(!is.finite(adj))) stop("non-finite heteroscedastic adjustment")
      adj
    }, error = function(e) NULL)
    if (!is.null(corrected)) {
      colnames(corrected) <- colnames(param_pick)
      return(corrected)
    }
  }

  adj <- sweep(residuals, 2, pred_at_target, "+", check.margin = FALSE)
  colnames(adj) <- colnames(param_pick)
  adj
}

fast_abc_estimate_one <- function(target, param, sumstat, tol, estimator = NULL,
                                  hcorr = TRUE) {
  estimator <- fast_estimator_mode(estimator)
  if (is.vector(sumstat)) sumstat <- matrix(sumstat, ncol = 1)
  if (is.vector(param)) param <- matrix(param, ncol = 1)

  scales <- mad_scales(sumstat)
  scaled_sumstat <- sweep(sumstat, 2, scales, "/", check.margin = FALSE)
  scaled_target <- target / scales
  d2 <- rowSums(sweep(scaled_sumstat, 2, scaled_target, "-", check.margin = FALSE)^2)
  pick <- accepted_indices_sq(d2, tol)
  if (length(pick) < 2L) return(NA_real_)

  param_pick <- param[pick, , drop = FALSE]
  if (estimator == "rejection") return(qst_from_param_matrix(param_pick))

  weights <- epanechnikov_weights_sq(d2[pick])
  adj <- tryCatch(
    loclinear_adjust(scaled_target, param_pick, scaled_sumstat[pick, , drop = FALSE],
                     weights, hcorr = hcorr),
    error = function(e) NULL
  )
  if (is.null(adj)) return(qst_from_param_matrix(param_pick))
  qst_from_param_matrix(adj)
}

combo_name <- function(combo) paste(combo, collapse = ",")

build_combo_indices <- function(summary_stat_combos, stat_names) {
  idx_list <- lapply(summary_stat_combos, function(combo) match(combo, stat_names))
  if (any(vapply(idx_list, function(x) any(is.na(x)), logical(1)))) {
    stop("At least one combo contains a stat missing from simulated stats")
  }
  idx_list
}

fill_chunk_mask <- function(combo_indices_chunk, n_cols) {
  rows <- rep.int(seq_along(combo_indices_chunk), lengths(combo_indices_chunk))
  cols <- unlist(combo_indices_chunk, use.names = FALSE)
  mask <- matrix(0, nrow = length(combo_indices_chunk), ncol = n_cols)
  mask[cbind(rows, cols)] <- 1
  mask
}

fast_abc_estimate_multi <- function(obs_stats, prior_params_valid, sim_stats_valid,
                                    summary_stat_combos, tol, estimator = NULL,
                                    chunk_size = NULL, hcorr = TRUE) {
  estimator <- fast_estimator_mode(estimator)
  if (is.null(chunk_size)) {
    chunk_size <- as.integer(Sys.getenv("FAST_ABC_COMBO_CHUNK_SIZE", "64"))
  }
  chunk_size <- max(1L, chunk_size)

  stat_names <- colnames(sim_stats_valid)
  scales <- mad_scales(sim_stats_valid)
  scaled_sim <- sweep(sim_stats_valid, 2, scales, "/", check.margin = FALSE)
  scaled_target <- obs_stats[stat_names] / scales
  diff2 <- sweep(scaled_sim, 2, scaled_target, "-", check.margin = FALSE)^2

  combo_indices <- build_combo_indices(summary_stat_combos, stat_names)
  combo_names <- vapply(summary_stat_combos, combo_name, character(1))
  total <- length(summary_stat_combos)
  out <- rep(NA_real_, total)
  names(out) <- combo_names

  sb <- prior_params_valid[, "sd_between_pop"]
  sw <- prior_params_valid[, "sd_within_pop"]
  sb2 <- sb * sb
  qst_pool <- sb2 / (sb2 + 2 * sw * sw)

  for (start in seq(1L, total, by = chunk_size)) {
    end <- min(total, start + chunk_size - 1L)
    mask <- fill_chunk_mask(combo_indices[start:end], ncol(diff2))
    dist2_chunk <- tcrossprod(diff2, mask)
    chunk_len <- end - start + 1L

    for (j in seq_len(chunk_len)) {
      i <- start + j - 1L
      d2 <- dist2_chunk[, j]
      pick <- accepted_indices_sq(d2, tol)
      if (length(pick) < 2L) next

      if (estimator == "rejection") {
        out[i] <- mean(qst_pool[pick])
        next
      }

      idx <- combo_indices[[i]]
      param_pick <- prior_params_valid[pick, , drop = FALSE]
      weights <- epanechnikov_weights_sq(d2[pick])
      adj <- tryCatch(
        loclinear_adjust(
          scaled_target[idx],
          param_pick,
          scaled_sim[pick, idx, drop = FALSE],
          weights,
          hcorr = hcorr
        ),
        error = function(e) NULL
      )
      out[i] <- if (is.null(adj)) qst_from_param_matrix(param_pick) else qst_from_param_matrix(adj)
    }

    rm(dist2_chunk, mask)
  }

  gc(verbose = FALSE, reset = TRUE)
  out
}

build_fast_reference_pool <- function(obs_stats, num_sim, summary_stat_combos) {
  required_stats <- get_required_stats(summary_stat_combos)
  n_stats <- length(required_stats)
  batch_size <- 5000L
  num_batches <- ceiling(num_sim / batch_size)

  prior_floor <- 0.1 * (obs_stats["among_pop_sd"] + obs_stats["within_pop_sd"] + obs_stats["ext_sd"])
  sd_between_pop_prior <- runif(num_sim, 0, max(2 * obs_stats["among_pop_sd"], prior_floor))
  sd_within_pop_prior <- runif(num_sim, 0, max(2 * obs_stats["within_pop_sd"], prior_floor))
  sd_ext_prior <- runif(num_sim, 0, 2 * obs_stats["ext_sd"])

  sim_stats_matrix <- matrix(NA_real_, nrow = num_sim, ncol = n_stats)
  prior_params <- cbind(sd_between_pop_prior, sd_within_pop_prior, sd_ext_prior)
  colnames(prior_params) <- c("sd_between_pop", "sd_within_pop", "sd_ext")
  rm(sd_between_pop_prior, sd_within_pop_prior, sd_ext_prior)
  gc(verbose = FALSE, reset = TRUE)

  for (b in seq_len(num_batches)) {
    start_idx <- (b - 1L) * batch_size + 1L
    end_idx <- min(b * batch_size, num_sim)
    batch_results <- run_batch_simulations(
      end_idx - start_idx + 1L,
      prior_params[start_idx:end_idx, 1],
      prior_params[start_idx:end_idx, 2],
      prior_params[start_idx:end_idx, 3],
      required_stats = required_stats
    )
    sim_stats_matrix[start_idx:end_idx, ] <- batch_results
    rm(batch_results)
    gc(verbose = FALSE, reset = TRUE)
  }

  colnames(sim_stats_matrix) <- required_stats
  valid_rows <- complete.cases(sim_stats_matrix) & (rowSums(is.infinite(sim_stats_matrix)) == 0)
  sim_stats_valid <- sim_stats_matrix[valid_rows, , drop = FALSE]
  prior_params_valid <- prior_params[valid_rows, , drop = FALSE]
  rm(sim_stats_matrix, prior_params, valid_rows)
  gc(verbose = FALSE, reset = TRUE)

  list(
    sim_stats_valid = sim_stats_valid,
    prior_params_valid = prior_params_valid,
    required_stats = required_stats
  )
}

estimate_qst_abc_all_combos_fast <- function(obs_stats, num_sim, summary_stat_combos,
                                             estimator = NULL) {
  pool <- build_fast_reference_pool(obs_stats, num_sim, summary_stat_combos)
  n_valid <- nrow(pool$sim_stats_valid)
  tol <- max(0.001, get_abc_tol_numerator() / n_valid)
  cat("Fast ABC estimator (v2):", fast_estimator_mode(estimator), "\n")
  cat("Valid simulations:", n_valid, "/", num_sim, "\n")
  fast_abc_estimate_multi(
    obs_stats = obs_stats[pool$required_stats],
    prior_params_valid = pool$prior_params_valid,
    sim_stats_valid = pool$sim_stats_valid,
    summary_stat_combos = summary_stat_combos,
    tol = tol,
    estimator = estimator
  )
}

process_evaluate_batch_all_combos_fast <- function(n_repeats, qst_value, ve_ratio, num_sim,
                                                   summary_stat_combos, start_id = 1L,
                                                   estimator = NULL) {
  n_combos <- length(summary_stat_combos)
  combo_names <- vapply(summary_stat_combos, combo_name, character(1))
  required_stats <- get_required_stats(summary_stat_combos)
  results <- data.frame(
    repeat_id = rep(seq(start_id, start_id + n_repeats - 1L), each = n_combos),
    true_qst = rep(qst_value, n_repeats * n_combos),
    ve_ratio = rep(ve_ratio, n_repeats * n_combos),
    combo = rep(combo_names, n_repeats),
    estimated_qst = rep(NA_real_, n_repeats * n_combos),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_repeats)) {
    repeat_id <- start_id + i - 1L
    seed <- (round(qst_value * 100) * 10000 + round(ve_ratio * 1000) * 100 + repeat_id) %%
      .Machine$integer.max
    obs_stats <- generate_adaptive_obs_stats(qst_value, ve_ratio, seed, required_stats)
    if (!is.null(obs_stats)) {
      qst_estimates <- estimate_qst_abc_all_combos_fast(
        obs_stats, num_sim, summary_stat_combos, estimator = estimator
      )
      base <- (i - 1L) * n_combos
      results$estimated_qst[base + seq_len(n_combos)] <- as.numeric(qst_estimates[combo_names])
    }
    gc(verbose = FALSE, full = TRUE, reset = TRUE)
    if (i %% 10 == 0 || i == n_repeats) {
      cat("  Processed", i, "/", n_repeats, "repeats\n")
    }
  }
  results
}

process_fst_batch_all_combos_fast <- function(fst_values, ratioVext, num_sim,
                                              summary_stat_combos, start_idx = 1L,
                                              estimator = NULL) {
  n_fst <- length(fst_values)
  n_combos <- length(summary_stat_combos)
  combo_names <- vapply(summary_stat_combos, combo_name, character(1))
  required_stats <- get_required_stats(summary_stat_combos)
  results <- data.frame(
    fst = rep(fst_values, each = n_combos),
    combo = rep(combo_names, n_fst),
    qst = rep(NA_real_, n_fst * n_combos),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(fst_values)) {
    fst_value <- fst_values[i]
    global_idx <- start_idx + i - 1L
    seed <- (round(fst_value * 1e6) + global_idx * 10000) %% .Machine$integer.max
    obs_stats <- generate_neutral_obs_stats(fst_value, ratioVext, seed, required_stats)
    if (!is.null(obs_stats)) {
      qst_estimates <- estimate_qst_abc_all_combos_fast(
        obs_stats, num_sim, summary_stat_combos, estimator = estimator
      )
      base <- (i - 1L) * n_combos
      results$qst[base + seq_len(n_combos)] <- as.numeric(qst_estimates[combo_names])
    }
    gc(verbose = FALSE, full = TRUE, reset = TRUE)
    if (i %% 10 == 0 || i == n_fst) {
      cat("  Processed", i, "/", n_fst, "FST values\n")
    }
  }
  results
}

run_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 6) {
    cat("Usage:\n")
    cat("  adaptive: Rscript perf_eval_all_combos_fast_cell_v2.R adaptive START N QST VE combos out [num_sim]\n")
    cat("  neutral:  Rscript perf_eval_all_combos_fast_cell_v2.R neutral FST_FILE START_IDX VE combos out [num_sim]\n")
    quit(status = 1)
  }

  mode <- args[[1]]
  num_sim_default <- 100000L
  estimator <- fast_estimator_mode()

  if (mode == "adaptive") {
    if (length(args) < 7) stop("adaptive needs 7+ args")
    start_id <- as.integer(args[[2]])
    n_repeats <- as.integer(args[[3]])
    qst_value <- as.numeric(args[[4]])
    ve_ratio <- as.numeric(args[[5]])
    combos_file <- args[[6]]
    output_file <- args[[7]]
    num_sim <- if (length(args) >= 8) as.integer(args[[8]]) else num_sim_default
    summary_stat_combos <- parse_summary_stats(combos_file)
    cat("Fast all-combo adaptive cell (v2): QST=", qst_value, " VE=", ve_ratio,
        " repeats=", n_repeats, " start=", start_id, " num_sim=", num_sim,
        " combos=", length(summary_stat_combos), " estimator=", estimator, "\n", sep = "")
    batch_results <- process_evaluate_batch_all_combos_fast(
      n_repeats, qst_value, ve_ratio, num_sim, summary_stat_combos, start_id, estimator
    )
    result <- list(
      mode = "batch_evaluate",
      true_qst = qst_value,
      ve_ratio = ve_ratio,
      n_repeats = n_repeats,
      start_id = start_id,
      is_multi_combo = TRUE,
      results = batch_results,
      fast_all_combos = list(num_sim = num_sim, combos_file = combos_file,
                             estimator = estimator, version = "v2")
    )
  } else if (mode == "neutral") {
    if (length(args) < 7) stop("neutral needs 7+ args")
    fst_file <- args[[2]]
    start_idx <- as.integer(args[[3]])
    ve_ratio <- as.numeric(args[[4]])
    combos_file <- args[[5]]
    output_file <- args[[6]]
    num_sim <- if (length(args) >= 7) as.integer(args[[7]]) else num_sim_default
    fst_values <- scan(fst_file, what = numeric(), quiet = TRUE)
    summary_stat_combos <- parse_summary_stats(combos_file)
    cat("Fast all-combo neutral cell (v2): n_fst=", length(fst_values), " start_idx=", start_idx,
        " VE=", ve_ratio, " num_sim=", num_sim, " combos=", length(summary_stat_combos),
        " estimator=", estimator, "\n", sep = "")
    batch_results <- process_fst_batch_all_combos_fast(
      fst_values, ve_ratio, num_sim, summary_stat_combos, start_idx, estimator
    )
    result <- list(
      mode = "batch_neutral",
      n_fst = length(fst_values),
      is_multi_combo = TRUE,
      results = batch_results,
      fast_all_combos = list(num_sim = num_sim, combos_file = combos_file,
                             estimator = estimator, version = "v2")
    )
  } else {
    stop("mode must be adaptive or neutral, got: ", mode)
  }

  save(result, file = output_file)
  cat("Saved:", output_file, "\n")
}

if (!interactive() && Sys.getenv("PERF_EVAL_ALL_COMBOS_FAST_SOURCE_ONLY", "") != "1") {
  run_cli()
}
