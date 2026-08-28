# Vectorized local-linear and rejection ABC for QST (matches abc::abc loclinear
# on the same reference pool; used when QST_ABC_METHOD is loclinear or rejection).
# Sourced from qst_abc_sim.R (same directory). Neural network / ridge still use abc::abc.

resolve_fast_estimator_arg <- function(estimator) {
  if (is.null(estimator)) {
    estimator <- trimws(Sys.getenv("FAST_ABC_ESTIMATOR", "loclinear"))
  }
  if (!estimator %in% c("loclinear", "rejection")) {
    stop("fast ABC estimator must be loclinear or rejection, got: ", estimator)
  }
  estimator
}

mad_scale <- function(x) {
  s <- stats::mad(x)
  if (!is.finite(s) || s == 0) 1 else s
}

mad_scales <- function(mat) {
  vapply(seq_len(ncol(mat)), function(j) mad_scale(mat[, j]), numeric(1))
}

accepted_indices_sq <- function(d2, tol) {
  n <- length(d2)
  nacc <- max(1L, min(n, ceiling(n * tol)))
  ds <- sort(d2, partial = nacc)[nacc]
  idx <- which(d2 <= ds)
  if (length(idx) > nacc) {
    idx <- idx[1:nacc]
  }
  idx
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

qst_per_particle <- function(param) {
  sb <- param[, "sd_between_pop"]
  sw <- param[, "sd_within_pop"]
  sb2 <- sb * sb
  sb2 / (sb2 + 2 * sw * sw)
}

#' Chromosome F_ST 95th-percentile reference for D1 tail detector (Design Module).
get_qst_tail_threshold <- function() {
  v <- suppressWarnings(as.numeric(Sys.getenv("QST_TAIL_THRESHOLD_FST95", "")))
  if (length(v) == 1L && is.finite(v)) return(v)
  chr <- trimws(Sys.getenv("QST_PERF_EVAL_CHR", "autosomes"))
  if (identical(chr, "chrX")) return(0.705882)
  0.615825
}

p_qst_gt_threshold_from_param <- function(param, threshold = NULL) {
  if (is.null(param) || nrow(param) < 1L) return(NA_real_)
  if (is.null(threshold)) threshold <- get_qst_tail_threshold()
  mean(qst_per_particle(param) > threshold, na.rm = TRUE)
}

p_qst_gt_threshold_from_qst_vec <- function(qst_vec, threshold = NULL) {
  if (is.null(qst_vec) || !length(qst_vec)) return(NA_real_)
  if (is.null(threshold)) threshold <- get_qst_tail_threshold()
  mean(qst_vec > threshold, na.rm = TRUE)
}

#' Extract QST point estimates from fast_abc_estimate_multi return value.
fast_abc_multi_qst <- function(res) {
  if (is.list(res) && !is.null(res$qst)) return(res$qst)
  res
}

fast_abc_multi_p_tail <- function(res) {
  if (is.list(res) && !is.null(res$p_qst_gt_fst95)) return(res$p_qst_gt_fst95)
  rep(NA_real_, length(fast_abc_multi_qst(res)))
}

ratioVext_per_particle <- function(param) {
  vgb <- param[, "sd_between_pop"]^2
  vgw <- param[, "sd_within_pop"]^2
  ve <- param[, "sd_ext"]^2
  vg <- (vgb + 2 * vgw) / 2
  ifelse(is.finite(vg) & vg > 0, ve / vg, NA_real_)
}

abc_posterior_quantiles_from_param <- function(param, probs = c(0.025, 0.975)) {
  if (is.null(param) || nrow(param) < 1L) {
    return(list(
      qst_lo = NA_real_, qst_hi = NA_real_,
      ratioVext_lo = NA_real_, ratioVext_hi = NA_real_
    ))
  }
  q_qst <- stats::quantile(qst_per_particle(param), probs = probs, na.rm = TRUE, names = FALSE)
  q_rv <- stats::quantile(ratioVext_per_particle(param), probs = probs, na.rm = TRUE, names = FALSE)
  list(
    qst_lo = as.numeric(q_qst[1]),
    qst_hi = as.numeric(q_qst[2]),
    ratioVext_lo = as.numeric(q_rv[1]),
    ratioVext_hi = as.numeric(q_rv[2])
  )
}

fast_abc_result_from_param <- function(param, threshold = NULL) {
  q <- abc_posterior_quantiles_from_param(param)
  sds <- abc_sd_means_from_param(param)
  c(
    list(
      qst = qst_from_param_matrix(param),
      p_qst_gt_fst95 = p_qst_gt_threshold_from_param(param, threshold),
      qst_lo = q$qst_lo,
      qst_hi = q$qst_hi,
      ratioVext_lo = q$ratioVext_lo,
      ratioVext_hi = q$ratioVext_hi
    ),
    sds
  )
}

abc_sd_means_from_param <- function(param) {
  if (is.null(param) || nrow(param) < 1L) {
    return(list(
      abc_among_pop_sd = NA_real_,
      abc_within_pop_sd = NA_real_,
      abc_ext_sd = NA_real_
    ))
  }
  list(
    abc_among_pop_sd = mean(param[, "sd_between_pop"], na.rm = TRUE),
    abc_within_pop_sd = mean(param[, "sd_within_pop"], na.rm = TRUE),
    abc_ext_sd = mean(param[, "sd_ext"], na.rm = TRUE)
  )
}

abc_sd_means_from_adj <- function(adj) {
  abc_sd_means_from_param(adj)
}

abc_ratioVext_from_sds <- function(among_sd, within_sd, ext_sd) {
  if (!all(is.finite(c(among_sd, within_sd, ext_sd)))) return(NA_real_)
  vg <- (among_sd^2 + 2 * within_sd^2) / 2
  if (!is.finite(vg) || vg <= 0) return(NA_real_)
  ext_sd^2 / vg
}

calib_extended_output_enabled <- function() {
  identical(toupper(Sys.getenv("QST_CALIB_EXTENDED_OUTPUT", "FALSE")), "TRUE")
}

fst_batch_row_from_abc <- function(fst_value, res, both_var_negative = 0L) {
  row <- list(
    fst = fst_value,
    qst = res$qst,
    both_var_negative = as.integer(both_var_negative)
  )
  if (!calib_extended_output_enabled()) return(row)
  among <- as.numeric(res$abc_among_pop_sd)
  within <- as.numeric(res$abc_within_pop_sd)
  ext <- as.numeric(res$abc_ext_sd)
  c(
    row,
    list(
      abc_among_pop_sd = among,
      abc_within_pop_sd = within,
      abc_ext_sd = ext,
      abc_ratioVext = abc_ratioVext_from_sds(among, within, ext),
      qst_lo = as.numeric(res$qst_lo),
      qst_hi = as.numeric(res$qst_hi)
    )
  )
}

fast_abc_empty_result <- function() {
  c(
    list(
      qst = NA_real_,
      p_qst_gt_fst95 = NA_real_,
      qst_lo = NA_real_,
      qst_hi = NA_real_,
      ratioVext_lo = NA_real_,
      ratioVext_hi = NA_real_
    ),
    abc_sd_means_from_param(NULL)
  )
}

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
  for (j in seq_len(ncol(residuals))) {
    residuals[, j] <- residuals[, j] - resid_mean[j]
  }
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
      for (j in seq_len(ncol(tmp))) {
        tmp[, j] <- tmp[, j] * pred_sd_vec[j]
      }
      adj <- tmp
      for (j in seq_len(ncol(adj))) {
        adj[, j] <- adj[, j] + pred_at_target[j]
      }
      if (any(!is.finite(adj))) stop("non-finite heteroscedastic adjustment")
      adj
    }, error = function(e) NULL)
    if (!is.null(corrected)) {
      colnames(corrected) <- colnames(param_pick)
      return(corrected)
    }
  }

  adj <- residuals
  for (j in seq_len(ncol(adj))) {
    adj[, j] <- adj[, j] + pred_at_target[j]
  }
  colnames(adj) <- colnames(param_pick)
  adj
}

fast_abc_estimate_one <- function(target, param, sumstat, tol, estimator = NULL,
                                  hcorr = TRUE) {
  estimator <- resolve_fast_estimator_arg(estimator)
  if (is.vector(sumstat)) sumstat <- matrix(sumstat, ncol = 1)
  if (is.vector(param)) param <- matrix(param, ncol = 1)

  scales <- mad_scales(sumstat)
  scaled_sumstat <- sumstat
  for (j in seq_len(ncol(sumstat))) {
    scaled_sumstat[, j] <- sumstat[, j] / scales[j]
  }
  scaled_target <- target / scales
  d2 <- 0
  for (j in seq_len(ncol(scaled_sumstat))) {
    d2 <- d2 + (scaled_sumstat[, j] - scaled_target[j])^2
  }
  pick <- accepted_indices_sq(d2, tol)
  if (length(pick) < 2L) return(fast_abc_empty_result())

  param_pick <- param[pick, , drop = FALSE]
  if (estimator == "rejection") {
    return(fast_abc_result_from_param(param_pick))
  }

  weights <- epanechnikov_weights_sq(d2[pick])
  adj <- tryCatch(
    loclinear_adjust(scaled_target, param_pick, scaled_sumstat[pick, , drop = FALSE],
                     weights, hcorr = hcorr),
    error = function(e) NULL
  )
  if (is.null(adj)) {
    return(fast_abc_result_from_param(param_pick))
  }
  fast_abc_result_from_param(adj)
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
  estimator <- resolve_fast_estimator_arg(estimator)
  if (is.null(chunk_size)) {
    chunk_size <- as.integer(Sys.getenv("FAST_ABC_COMBO_CHUNK_SIZE", "64"))
  }
  chunk_size <- max(1L, chunk_size)

  stat_names <- colnames(sim_stats_valid)
  scales <- mad_scales(sim_stats_valid)
  scaled_sim <- sim_stats_valid
  for (j in seq_len(ncol(sim_stats_valid))) {
    scaled_sim[, j] <- sim_stats_valid[, j] / scales[j]
  }
  scaled_target <- obs_stats[stat_names] / scales
  diff2 <- scaled_sim
  for (j in seq_len(ncol(scaled_sim))) {
    diff2[, j] <- (scaled_sim[, j] - scaled_target[j])^2
  }

  combo_indices <- build_combo_indices(summary_stat_combos, stat_names)
  combo_names <- vapply(summary_stat_combos, combo_name, character(1))
  total <- length(summary_stat_combos)
  out <- rep(NA_real_, total)
  out_p <- rep(NA_real_, total)
  names(out) <- combo_names
  names(out_p) <- combo_names
  tail_thr <- get_qst_tail_threshold()

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
        out_p[i] <- p_qst_gt_threshold_from_qst_vec(qst_pool[pick], tail_thr)
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
      param_use <- if (is.null(adj)) param_pick else adj
      out[i] <- qst_from_param_matrix(param_use)
      out_p[i] <- p_qst_gt_threshold_from_param(param_use, tail_thr)
    }

    rm(dist2_chunk, mask)
  }

  gc(verbose = FALSE, reset = TRUE)
  list(qst = out, p_qst_gt_fst95 = out_p)
}

fast_process_fst_batch <- function(fst_values, ratioVext, num_sim, summary_stat_names = c("QST", "ratioVext"), start_idx = 1) {
  n_fst <- length(fst_values)
  
  condor_cores <- suppressWarnings(as.integer(Sys.getenv("_CONDOR_NPROCS", "1")))
  num_cores <- suppressWarnings(as.integer(Sys.getenv("QST_ABC_CORES", as.character(condor_cores))))
  if (is.na(num_cores) || num_cores < 1L) num_cores <- 1L
  num_cores <- min(num_cores, parallel::detectCores())
  
  cat("Running FST batch processing in parallel using", num_cores, "cores\n")
  
  results_list <- parallel::mclapply(seq_along(fst_values), function(i) {
    fst_value <- fst_values[i]
    global_idx <- start_idx + i - 1
    seed <- (round(fst_value * 1e6) + global_idx * 10000) %% .Machine$integer.max
    obs <- generate_neutral_obs_stats(fst_value, ratioVext, seed, c(summary_stat_names, "among_pop_sd", "within_pop_sd", "ext_sd"))
    
    if (is.null(obs)) {
      return(fst_batch_row_from_abc(fst_value, fast_abc_empty_result(), 0L))
    }
    if (is_both_neg(obs)) {
      return(fst_batch_row_from_abc(fst_value, fast_abc_empty_result(), 1L))
    }
    
    prior_floor <- 0.1 * (obs['among_pop_sd'] + obs['within_pop_sd'])
    max_sd_between <- max(2 * obs["among_pop_sd"], prior_floor)
    max_sd_within <- max(2 * obs["within_pop_sd"], prior_floor)
    max_sd_ext <- 2 * obs["ext_sd"]
    
    abc_seed <- abc_seed_from_obs_stats(obs)
    set.seed(abc_seed)
    
    sd_between_pop_prior <- runif(num_sim, 0, max_sd_between)
    sd_within_pop_prior  <- runif(num_sim, 0, max_sd_within)
    sd_ext_prior         <- runif(num_sim, 0, max_sd_ext)
    prior_params <- cbind(sd_between_pop_prior, sd_within_pop_prior, sd_ext_prior)
    colnames(prior_params) <- c("sd_between_pop", "sd_within_pop", "sd_ext")
    rm(sd_between_pop_prior, sd_within_pop_prior, sd_ext_prior)
    
    batch_size <- as.integer(Sys.getenv("QST_ABC_BATCH_SIZE", "50000"))
    num_batches <- ceiling(num_sim / batch_size)
    sim_stats_matrix <- matrix(NA, nrow = num_sim, ncol = length(summary_stat_names))
    colnames(sim_stats_matrix) <- summary_stat_names
    
    for (b in 1:num_batches) {
      s_idx <- (b - 1) * batch_size + 1
      e_idx <- min(b * batch_size, num_sim)
      batch_results <- run_batch_simulations(
        e_idx - s_idx + 1,
        prior_params[s_idx:e_idx, 1],
        prior_params[s_idx:e_idx, 2],
        prior_params[s_idx:e_idx, 3],
        required_stats = summary_stat_names
      )
      sim_stats_matrix[s_idx:e_idx, ] <- batch_results
      rm(batch_results)
    }
    
    valid_rows <- complete.cases(sim_stats_matrix) & (rowSums(is.infinite(sim_stats_matrix)) == 0)
    
    if (sum(valid_rows) < 2L) {
      return(fst_batch_row_from_abc(fst_value, fast_abc_empty_result(), 0L))
    }
    
    sim_stats_valid <- sim_stats_matrix[valid_rows, , drop = FALSE]
    prior_params_valid <- prior_params[valid_rows, , drop = FALSE]
    
    res <- fast_abc_estimate_one(
      obs[summary_stat_names],
      prior_params_valid,
      sim_stats_valid,
      get_abc_tol_numerator() / nrow(sim_stats_valid),
      get_abc_method()
    )
    return(fst_batch_row_from_abc(fst_value, res, 0L))
  }, mc.cores = num_cores)
  
  results <- do.call(rbind, lapply(results_list, function(x) {
    as.data.frame(x, stringsAsFactors = FALSE)
  }))
  return(results)
}

fast_process_fst_batch_multi <- function(fst_values, ratioVext, num_sim, summary_stat_combos, start_idx = 1) {
  n_fst <- length(fst_values)
  n_combos <- length(summary_stat_combos)
  combo_names <- sapply(summary_stat_combos, function(x) paste(x, collapse = ","))
  
  condor_cores <- suppressWarnings(as.integer(Sys.getenv("_CONDOR_NPROCS", "1")))
  num_cores <- suppressWarnings(as.integer(Sys.getenv("QST_ABC_CORES", as.character(condor_cores))))
  if (is.na(num_cores) || num_cores < 1L) num_cores <- 1L
  num_cores <- min(num_cores, parallel::detectCores())
  
  cat("Running FST batch processing multi in parallel using", num_cores, "cores\n")
  
  all_required_stats <- unique(unlist(summary_stat_combos))
  
  results_list <- parallel::mclapply(seq_along(fst_values), function(i) {
    fst_value <- fst_values[i]
    global_idx <- start_idx + i - 1
    seed <- (round(fst_value * 1e6) + global_idx * 10000) %% .Machine$integer.max
    obs <- generate_neutral_obs_stats(fst_value, ratioVext, seed, unique(c(all_required_stats, "among_pop_sd", "within_pop_sd", "ext_sd")))
    
    res_df <- data.frame(
      fst = rep(fst_value, n_combos),
      combo = combo_names,
      qst = rep(NA_real_, n_combos),
      p_qst_gt_fst95 = rep(NA_real_, n_combos),
      both_var_negative = rep(0L, n_combos),
      stringsAsFactors = FALSE
    )
    
    if (is.null(obs)) {
      return(res_df)
    }
    if (is_both_neg(obs)) {
      res_df$both_var_negative <- 1L
      return(res_df)
    }
    
    prior_floor <- 0.1 * (obs['among_pop_sd'] + obs['within_pop_sd'])
    max_sd_between <- max(2 * obs["among_pop_sd"], prior_floor)
    max_sd_within <- max(2 * obs["within_pop_sd"], prior_floor)
    max_sd_ext <- 2 * obs["ext_sd"]
    
    abc_seed <- abc_seed_from_obs_stats(obs)
    set.seed(abc_seed)
    
    sd_between_pop_prior <- runif(num_sim, 0, max_sd_between)
    sd_within_pop_prior  <- runif(num_sim, 0, max_sd_within)
    sd_ext_prior         <- runif(num_sim, 0, max_sd_ext)
    prior_params <- cbind(sd_between_pop_prior, sd_within_pop_prior, sd_ext_prior)
    colnames(prior_params) <- c("sd_between_pop", "sd_within_pop", "sd_ext")
    rm(sd_between_pop_prior, sd_within_pop_prior, sd_ext_prior)
    
    batch_size <- as.integer(Sys.getenv("QST_ABC_BATCH_SIZE", "50000"))
    num_batches <- ceiling(num_sim / batch_size)
    sim_stats_matrix <- matrix(NA, nrow = num_sim, ncol = length(all_required_stats))
    colnames(sim_stats_matrix) <- all_required_stats
    
    for (b in 1:num_batches) {
      s_idx <- (b - 1) * batch_size + 1
      e_idx <- min(b * batch_size, num_sim)
      batch_results <- run_batch_simulations(
        e_idx - s_idx + 1,
        prior_params[s_idx:e_idx, 1],
        prior_params[s_idx:e_idx, 2],
        prior_params[s_idx:e_idx, 3],
        required_stats = all_required_stats
      )
      sim_stats_matrix[s_idx:e_idx, ] <- batch_results
      rm(batch_results)
    }
    
    valid_rows <- complete.cases(sim_stats_matrix) & (rowSums(is.infinite(sim_stats_matrix)) == 0)
    
    if (sum(valid_rows) < 2L) {
      return(res_df)
    }
    
    sim_stats_valid <- sim_stats_matrix[valid_rows, , drop = FALSE]
    prior_params_valid <- prior_params[valid_rows, , drop = FALSE]
    
    res <- fast_abc_estimate_multi(
      obs,
      prior_params_valid,
      sim_stats_valid,
      summary_stat_combos,
      get_abc_tol_numerator() / nrow(sim_stats_valid),
      get_abc_method()
    )
    res_df$qst <- fast_abc_multi_qst(res)
    res_df$p_qst_gt_fst95 <- fast_abc_multi_p_tail(res)
    return(res_df)
  }, mc.cores = num_cores)
  
  results <- do.call(rbind, results_list)
  return(results)
}
