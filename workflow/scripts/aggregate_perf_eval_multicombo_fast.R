#!/usr/bin/env Rscript
# Fast multi-combo aggregation for perf-eval outputs.
#
# Usage:
#   Rscript aggregate_perf_eval_multicombo_fast.R <results_dir> <chr_type> <adaptive_qst> <ve_ratios> <threshold_percentile> <output_file> [combo_list_file]
#
# Detector stats (env PERF_EVAL_DETECTOR_STATS, default "qst,p_qst_gt_fst95"):
#   qst              - ABC posterior mean QST (legacy point-estimate detector)
#   p_qst_gt_fst95   - D1: P(QST > chr F_ST 95th percentile | data)
#
# Writes per-combo CSVs per detector:
#   tpr_fpr_matrix_<chr>_combo_XXXX.csv                    (detector=qst)
#   tpr_fpr_matrix_<chr>_detector_p_qst_gt_fst95_combo_XXXX.csv

suppressPackageStartupMessages(library(parallel))
num_cores <- min(8L, max(1L, parallel::detectCores() - 1L))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript aggregate_perf_eval_multicombo_fast.R <results_dir> <chr_type> <adaptive_qst> <ve_ratios> <threshold> <output_file> [combo_list_file]")
}

results_dir <- args[1]
chr_type <- args[2]
adaptive_qst <- as.numeric(strsplit(args[3], ",")[[1]])
ve_ratios <- as.numeric(strsplit(args[4], ",")[[1]])
threshold_percentile <- as.numeric(args[5])
output_file <- args[6]
combo_list_file <- if (length(args) >= 7) args[7] else ""

detector_stats <- trimws(strsplit(Sys.getenv("PERF_EVAL_DETECTOR_STATS", "qst,p_qst_gt_fst95"), ",")[[1]])
detector_stats <- detector_stats[nzchar(detector_stats)]

read_combo_list_file <- function(path) {
  if (path == "" || !file.exists(path)) return(character(0))
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) return(character(0))
  gsub("\t", ",", lines, fixed = TRUE)
}

extract_rows <- function(result, detector) {
  if (is.null(result) || !is.list(result)) return(NULL)
  if (result$mode == "batch_neutral") {
    if (!is.null(result$is_multi_combo) && result$is_multi_combo && "combo" %in% names(result$results)) {
      d <- result$results
      val <- if (detector == "qst") {
        d$qst
      } else if (detector == "anova_qst") {
        if ("qst" %in% names(d)) d$qst else rep(NA_real_, nrow(d))
      } else if (detector == "anova_norep_qst" && "anova_norep_qst" %in% names(d)) {
        d$anova_norep_qst
      } else if (detector == "p_qst_gt_fst95" && "p_qst_gt_fst95" %in% names(d)) {
        d$p_qst_gt_fst95
      } else if (detector == "p_model_adaptive" && "p_model_adaptive" %in% names(d)) {
        d$p_model_adaptive
      } else {
        rep(NA_real_, nrow(d))
      }
      return(data.frame(combo = d$combo, value = val, stringsAsFactors = FALSE))
    }
    val <- if (detector == "qst") result$results$qst else NA_real_
    return(data.frame(combo = "default", value = val, stringsAsFactors = FALSE))
  }
  if (result$mode == "batch_evaluate") {
    if (!is.null(result$is_multi_combo) && result$is_multi_combo && "combo" %in% names(result$results)) {
      d <- result$results
      val <- if (detector == "qst") {
        d$estimated_qst
      } else if (detector == "anova_qst") {
        if ("estimated_qst" %in% names(d)) d$estimated_qst else rep(NA_real_, nrow(d))
      } else if (detector == "anova_norep_qst" && "anova_norep_qst" %in% names(d)) {
        d$anova_norep_qst
      } else if (detector == "p_qst_gt_fst95" && "p_qst_gt_fst95" %in% names(d)) {
        d$p_qst_gt_fst95
      } else if (detector == "p_model_adaptive" && "p_model_adaptive" %in% names(d)) {
        d$p_model_adaptive
      } else {
        rep(NA_real_, nrow(d))
      }
      return(data.frame(combo = d$combo, value = val, stringsAsFactors = FALSE))
    }
    val <- if (detector == "qst") result$results$estimated_qst else NA_real_
    return(data.frame(combo = "default", value = val, stringsAsFactors = FALSE))
  }
  NULL
}

load_estimates_dir <- function(dir_path, detector) {
  files <- list.files(dir_path, pattern = "\\.RData$", full.names = TRUE)
  if (length(files) == 0) {
    return(data.frame(combo = character(0), value = numeric(0), stringsAsFactors = FALSE))
  }
  out <- mclapply(files, function(f) {
    tryCatch({
      env <- new.env()
      load(f, envir = env)
      extract_rows(env$result, detector)
    }, error = function(e) {
      warning("Error loading ", f, ": ", e$message)
      NULL
    })
  }, mc.cores = num_cores)
  out <- Filter(Negate(is.null), out)
  if (length(out) == 0) {
    return(data.frame(combo = character(0), value = numeric(0), stringsAsFactors = FALSE))
  }
  combined <- do.call(rbind, out)
  combined <- combined[is.finite(combined$value), , drop = FALSE]
  combined
}

get_unique_combos_from_dir <- function(dir_path) {
  files <- list.files(dir_path, pattern = "\\.RData$", full.names = TRUE)
  for (f in files) {
    tryCatch({
      env <- new.env()
      load(f, envir = env)
      res <- env$result
      if (!is.null(res$is_multi_combo) && res$is_multi_combo && "combo" %in% names(res$results)) {
        return(sort(unique(as.character(res$results$combo))))
      }
    }, error = function(e) invisible(NULL))
  }
  character(0)
}

output_file_for_detector <- function(base_output, detector) {
  if (detector == "qst") return(base_output)
  sub("\\.csv$", paste0("_detector_", detector, ".csv"), base_output)
}

aggregate_one_detector <- function(detector, lev, n_combos, n_qst, n_ratios, adaptive_qst,
                                 neutral_cache, adaptive_cache, output_file, chr_type,
                                 threshold_percentile, ve_ratios) {
  threshold_mat <- matrix(NA_real_, nrow = n_combos, ncol = n_ratios)
  fpr_mat <- matrix(NA_real_, nrow = n_combos, ncol = n_ratios)
  tpr_arr <- array(NA_real_, dim = c(n_qst, n_ratios, n_combos))
  rownames(threshold_mat) <- lev
  rownames(fpr_mat) <- lev

  for (ri in seq_len(n_ratios)) {
    d <- neutral_cache[[ri]]
    if (!nrow(d)) next
    dc <- as.character(d$combo)
    row_idx_full <- match(dc, lev)
    okq <- !is.na(row_idx_full)
    th <- rep(NA_real_, n_combos)
    if (any(okq)) {
      tp <- tapply(d$value[okq], row_idx_full[okq], quantile, probs = threshold_percentile, na.rm = TRUE)
      th[as.integer(names(tp))] <- as.numeric(tp)
    }
    threshold_mat[, ri] <- th
    row_idx <- row_idx_full
    ok <- !is.na(row_idx)
    if (!any(ok)) next
    thr_per_row <- threshold_mat[row_idx[ok], ri]
    tp <- tapply(d$value[ok] > thr_per_row, row_idx[ok], mean)
    fpr_vec <- rep(NA_real_, n_combos)
    fpr_vec[as.integer(names(tp))] <- as.numeric(tp)
    fpr_mat[, ri] <- fpr_vec
  }

  for (qi in seq_len(n_qst)) {
    for (ri in seq_len(n_ratios)) {
      ad <- adaptive_cache[[qi]][[ri]]
      if (!nrow(ad)) next
      ac <- as.character(ad$combo)
      row_idx <- match(ac, lev)
      ok <- !is.na(row_idx)
      if (!any(ok)) next
      thr_per_row <- threshold_mat[row_idx[ok], ri]
      tp <- tapply(ad$value[ok] > thr_per_row, row_idx[ok], mean)
      tvec <- rep(NA_real_, n_combos)
      tvec[as.integer(names(tp))] <- as.numeric(tp)
      tpr_arr[qi, ri, ] <- tvec
    }
  }

  ve_cols <- paste0("VEratio_", ve_ratios)
  qst_row_names <- paste0("QST_", adaptive_qst)
  out_file <- output_file_for_detector(output_file, detector)
  output_dir <- dirname(out_file)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  mean_tpr <- apply(tpr_arr, 3, function(x) mean(x, na.rm = TRUE))
  mean_fpr <- rowMeans(fpr_mat, na.rm = TRUE)
  n_stats <- vapply(lev, function(s) length(strsplit(s, ",", fixed = TRUE)[[1]]), integer(1))

  model_ranking <- data.frame(
    model = lev,
    model_id = sprintf("combo_%04d", seq_len(n_combos)),
    detector = detector,
    mean_tpr = mean_tpr,
    mean_fpr = mean_fpr,
    n_stats = n_stats,
    stringsAsFactors = FALSE
  )
  model_ranking <- model_ranking[order(-model_ranking$mean_tpr), , drop = FALSE]
  model_ranking$rank <- seq_len(nrow(model_ranking))
  ranking_file <- sub("\\.csv$", "_model_ranking.csv", out_file)
  write.csv(model_ranking, ranking_file, row.names = FALSE)
  cat("Detector", detector, "- ranking:", ranking_file, "\n")

  for (ci in seq_len(n_combos)) {
    combo <- lev[ci]
    combo_id <- sprintf("combo_%04d", ci)
    thr_row <- as.numeric(threshold_mat[ci, , drop = TRUE])
    fpr_row <- as.numeric(fpr_mat[ci, , drop = TRUE])
    tpr_block <- tpr_arr[, , ci, drop = TRUE]
    base <- data.frame(
      type = c("threshold", "FPR", qst_row_names),
      rbind(matrix(thr_row, nrow = 1), matrix(fpr_row, nrow = 1), tpr_block),
      stringsAsFactors = FALSE
    )
    colnames(base)[-1] <- ve_cols
    base$model <- combo
    base$model_id <- combo_id
    base$detector <- detector
    combo_out <- sub("\\.csv$", paste0("_", combo_id, ".csv"), out_file)
    write.csv(base, combo_out, row.names = FALSE)
    if (detector == "qst" && !grepl(",", combo, fixed = TRUE)) {
      legacy <- base[, c("type", ve_cols), drop = FALSE]
      legacy_out <- file.path(output_dir, paste0("tpr_fpr_matrix_", chr_type, "_", combo, ".csv"))
      write.csv(legacy, legacy_out, row.names = FALSE)
    }
  }
  invisible(model_ranking)
}

cat("FAST multi-combo aggregation\n")
cat("Results directory:", results_dir, "\n")
cat("Chromosome type:", chr_type, "\n")
cat("Detectors:", paste(detector_stats, collapse = ", "), "\n")

n_qst <- length(adaptive_qst)
n_ratios <- length(ve_ratios)

first_neutral_dir <- file.path(results_dir, "neutral_ratio_1")
combo_probe <- get_unique_combos_from_dir(first_neutral_dir)
if (length(combo_probe) == 0) {
  stop("Fast path only supports multi-combo runs (expected combo column in RData).")
}

if (combo_list_file == "") {
  auto <- file.path(dirname(results_dir), "combinations_all_nonempty_11stats.txt")
  if (file.exists(auto)) combo_list_file <- auto
}
combo_list <- read_combo_list_file(combo_list_file)
if (length(combo_list) == 0) combo_list <- sort(unique(combo_probe))
n_combos <- length(combo_list)
lev <- combo_list
cat("Combos:", n_combos, "\n")

for (detector in detector_stats) {
  cat("Loading data for detector:", detector, "\n")
  neutral_cache <- vector("list", n_ratios)
  adaptive_cache <- vector("list", n_qst)
  for (qi in seq_len(n_qst)) adaptive_cache[[qi]] <- vector("list", n_ratios)

  for (ri in seq_len(n_ratios)) {
    nd <- file.path(results_dir, paste0("neutral_ratio_", ri))
    neutral_cache[[ri]] <- if (dir.exists(nd)) load_estimates_dir(nd, detector) else {
      data.frame(combo = character(0), value = numeric(0), stringsAsFactors = FALSE)
    }
  }
  for (qi in seq_len(n_qst)) {
    qst_str <- gsub("\\.", "_", sprintf("%.2f", adaptive_qst[qi]))
    for (ri in seq_len(n_ratios)) {
      ad <- file.path(results_dir, paste0("adaptive_q", qst_str, "_r", ri))
      adaptive_cache[[qi]][[ri]] <- if (dir.exists(ad)) load_estimates_dir(ad, detector) else {
        data.frame(combo = character(0), value = numeric(0), stringsAsFactors = FALSE)
      }
    }
  }

  aggregate_one_detector(
    detector, lev, n_combos, n_qst, n_ratios, adaptive_qst,
    neutral_cache, adaptive_cache, output_file, chr_type,
    threshold_percentile, ve_ratios
  )
}

cat("Done (fast). Detectors:", paste(detector_stats, collapse = ", "), "\n")
