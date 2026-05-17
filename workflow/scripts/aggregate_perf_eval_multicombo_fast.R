#!/usr/bin/env Rscript
# Fast multi-combo aggregation for perf-eval outputs.
#
# Loads each neutral_ratio_* / adaptive_q*_* directory once, then computes
# thresholds / FPR / TPR for all combos via merge + tapply (no per-combo RData reload).
#
# Usage (same as aggregate_perf_eval.R):
#   Rscript aggregate_perf_eval_multicombo_fast.R <results_dir> <chr_type> <adaptive_qst> <ve_ratios> <threshold_percentile> <output_file> [combo_list_file]
#
# Writes per-combo CSVs compatible with generate_combined_ranking.R:
#   tpr_fpr_matrix_<chr>_combo_XXXX.csv  (includes model + model_id columns)
#
# For single-token combos (no comma), also writes legacy 5-column CSV:
#   tpr_fpr_matrix_<chr>_<combo>.csv  (type + VEratio_* only; matches older naming)

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

read_combo_list_file <- function(path) {
  if (path == "" || !file.exists(path)) return(character(0))
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) return(character(0))
  gsub("\t", ",", lines, fixed = TRUE)
}

load_qst_estimates_dir <- function(dir_path) {
  files <- list.files(dir_path, pattern = "\\.RData$", full.names = TRUE)
  if (length(files) == 0) {
    return(data.frame(combo = character(0), qst = numeric(0), stringsAsFactors = FALSE))
  }
  out <- list()
  for (f in files) {
    tryCatch(
      {
        load(f)
        if (result$mode == "batch_neutral") {
          if (!is.null(result$is_multi_combo) && result$is_multi_combo && "combo" %in% names(result$results)) {
            out[[length(out) + 1]] <- result$results[, c("combo", "qst")]
          } else {
            out[[length(out) + 1]] <- data.frame(combo = "default", qst = result$results$qst, stringsAsFactors = FALSE)
          }
        } else if (result$mode == "batch_evaluate") {
          if (!is.null(result$is_multi_combo) && result$is_multi_combo && "combo" %in% names(result$results)) {
            out[[length(out) + 1]] <- data.frame(
              combo = result$results$combo,
              qst = result$results$estimated_qst,
              stringsAsFactors = FALSE
            )
          } else {
            out[[length(out) + 1]] <- data.frame(combo = "default", qst = result$results$estimated_qst, stringsAsFactors = FALSE)
          }
        }
      },
      error = function(e) warning("Error loading ", f, ": ", e$message)
    )
  }
  if (length(out) == 0) {
    return(data.frame(combo = character(0), qst = numeric(0), stringsAsFactors = FALSE))
  }
  combined <- do.call(rbind, out)
  combined <- combined[!is.na(combined$qst), , drop = FALSE]
  combined
}

get_unique_combos_from_dir <- function(dir_path) {
  files <- list.files(dir_path, pattern = "\\.RData$", full.names = TRUE)
  combos <- character(0)
  for (f in files) {
    tryCatch(
      {
        load(f)
        if (!is.null(result$is_multi_combo) && result$is_multi_combo && "combo" %in% names(result$results)) {
          combos <- unique(c(combos, as.character(result$results$combo)))
        }
      },
      error = function(e) invisible(NULL)
    )
  }
  combos
}

cat("FAST multi-combo aggregation\n")
cat("Results directory:", results_dir, "\n")
cat("Chromosome type:", chr_type, "\n")

n_qst <- length(adaptive_qst)
n_ratios <- length(ve_ratios)

first_neutral_dir <- file.path(results_dir, "neutral_ratio_1")
combo_probe <- get_unique_combos_from_dir(first_neutral_dir)
is_multi_combo <- length(combo_probe) > 0
if (!is_multi_combo) {
  stop("Fast path only supports multi-combo runs (expected combo column in RData).")
}

if (combo_list_file == "") {
  auto <- file.path(dirname(results_dir), "combinations_all_nonempty_11stats.txt")
  if (file.exists(auto)) combo_list_file <- auto
}
combo_list <- read_combo_list_file(combo_list_file)
if (length(combo_list) == 0) {
  combo_list <- sort(unique(combo_probe))
}
n_combos <- length(combo_list)
cat("Combos:", n_combos, "\n")

lev <- combo_list

neutral_cache <- vector("list", n_ratios)
adaptive_cache <- vector("list", n_qst)
for (qi in seq_len(n_qst)) adaptive_cache[[qi]] <- vector("list", n_ratios)

cat("Loading neutral batches...\n")
for (ri in seq_len(n_ratios)) {
  nd <- file.path(results_dir, paste0("neutral_ratio_", ri))
  neutral_cache[[ri]] <- if (dir.exists(nd)) load_qst_estimates_dir(nd) else data.frame(combo = character(0), qst = numeric(0))
}

cat("Loading adaptive batches...\n")
for (qi in seq_len(n_qst)) {
  qst_str <- gsub("\\.", "_", sprintf("%.2f", adaptive_qst[qi]))
  for (ri in seq_len(n_ratios)) {
    ad <- file.path(results_dir, paste0("adaptive_q", qst_str, "_r", ri))
    adaptive_cache[[qi]][[ri]] <- if (dir.exists(ad)) load_qst_estimates_dir(ad) else data.frame(combo = character(0), qst = numeric(0))
  }
}

threshold_mat <- matrix(NA_real_, nrow = n_combos, ncol = n_ratios)
fpr_mat <- matrix(NA_real_, nrow = n_combos, ncol = n_ratios)
tpr_arr <- array(NA_real_, dim = c(n_qst, n_ratios, n_combos))

rownames(threshold_mat) <- lev
rownames(fpr_mat) <- lev

# Avoid merge(d, thr_df): on multi-million-row d/ad it can allocate huge temporaries
# and segfault in merge.data.frame; match + numeric threshold is equivalent here.
cat("Computing thresholds and FPR...\n")
for (ri in seq_len(n_ratios)) {
  d <- neutral_cache[[ri]]
  if (!nrow(d)) next
  dc <- as.character(d$combo)
  row_idx_full <- match(dc, lev)
  okq <- !is.na(row_idx_full)
  th <- rep(NA_real_, n_combos)
  if (any(okq)) {
    tp <- tapply(d$qst[okq], row_idx_full[okq], quantile, probs = threshold_percentile, na.rm = TRUE)
    th[as.integer(names(tp))] <- as.numeric(tp)
  }
  threshold_mat[, ri] <- th

  row_idx <- row_idx_full
  ok <- !is.na(row_idx)
  if (!any(ok)) next
  thr_per_row <- threshold_mat[row_idx[ok], ri]
  tp <- tapply(d$qst[ok] > thr_per_row, row_idx[ok], mean)
  fpr_vec <- rep(NA_real_, n_combos)
  fpr_vec[as.integer(names(tp))] <- as.numeric(tp)
  fpr_mat[, ri] <- fpr_vec
}

cat("Computing TPR...\n")
for (qi in seq_len(n_qst)) {
  for (ri in seq_len(n_ratios)) {
    ad <- adaptive_cache[[qi]][[ri]]
    if (!nrow(ad)) next
    ac <- as.character(ad$combo)
    row_idx <- match(ac, lev)
    ok <- !is.na(row_idx)
    if (!any(ok)) next
    thr_per_row <- threshold_mat[row_idx[ok], ri]
    tp <- tapply(ad$qst[ok] > thr_per_row, row_idx[ok], mean)
    tvec <- rep(NA_real_, n_combos)
    tvec[as.integer(names(tp))] <- as.numeric(tp)
    tpr_arr[qi, ri, ] <- tvec
  }
}

ve_cols <- paste0("VEratio_", ve_ratios)
qst_row_names <- paste0("QST_", adaptive_qst)

output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

mean_tpr <- apply(tpr_arr, 3, function(x) mean(x, na.rm = TRUE))
mean_fpr <- rowMeans(fpr_mat, na.rm = TRUE)
n_stats <- vapply(lev, function(s) length(strsplit(s, ",", fixed = TRUE)[[1]]), integer(1))

model_ranking <- data.frame(
  model = lev,
  model_id = sprintf("combo_%04d", seq_len(n_combos)),
  mean_tpr = mean_tpr,
  mean_fpr = mean_fpr,
  n_stats = n_stats,
  stringsAsFactors = FALSE
)
model_ranking <- model_ranking[order(-model_ranking$mean_tpr), , drop = FALSE]
model_ranking$rank <- seq_len(nrow(model_ranking))
ranking_file <- sub("\\.csv$", "_model_ranking.csv", output_file)
write.csv(model_ranking, ranking_file, row.names = FALSE)
cat("Model ranking saved to:", ranking_file, "\n")

cat("Writing per-combo CSVs...\n")
for (ci in seq_len(n_combos)) {
  combo <- lev[ci]
  combo_id <- sprintf("combo_%04d", ci)

  thr_row <- as.numeric(threshold_mat[ci, , drop = TRUE])
  fpr_row <- as.numeric(fpr_mat[ci, , drop = TRUE])
  tpr_block <- tpr_arr[, , ci, drop = TRUE]

  base <- data.frame(
    type = c("threshold", "FPR", qst_row_names),
    rbind(
      matrix(thr_row, nrow = 1),
      matrix(fpr_row, nrow = 1),
      tpr_block
    ),
    stringsAsFactors = FALSE
  )
  colnames(base)[-1] <- ve_cols

  base$model <- combo
  base$model_id <- combo_id

  combo_out <- sub("\\.csv$", paste0("_", combo_id, ".csv"), output_file)
  write.csv(base, combo_out, row.names = FALSE)

  if (!grepl(",", combo, fixed = TRUE)) {
    legacy <- base[, c("type", ve_cols), drop = FALSE]
    legacy_out <- file.path(output_dir, paste0("tpr_fpr_matrix_", chr_type, "_", combo, ".csv"))
    write.csv(legacy, legacy_out, row.names = FALSE)
  }
}

cat("Done (fast). Wrote", n_combos, "combo CSVs (+ legacy single-token files).\n")
