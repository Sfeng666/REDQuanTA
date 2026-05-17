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
  estimator <- resolve_fast_estimator_arg(estimator)
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
  estimator <- resolve_fast_estimator_arg(estimator)
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
