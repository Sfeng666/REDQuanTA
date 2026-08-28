#!/usr/bin/env Rscript
# Aggregate neutral QST results and compare to trait QST

# Explicit both_var_negative flag only (no among==within soft-clip fallback).
is_both_neg <- function(obs_stats) {
  if ("both_var_negative" %in% names(obs_stats)) {
    v <- suppressWarnings(as.numeric(obs_stats["both_var_negative"]))
    if (length(v) == 1L && !is.na(v) && v == 1) return(TRUE)
  }
  FALSE
}

obs_num <- function(stats, key) {
  if (is.null(stats) || length(stats) == 0) return(NA_real_)
  v <- suppressWarnings(as.numeric(unlist(stats[key]))[1])
  if (length(v) == 0 || is.na(v)) NA_real_ else v
}

aggregate_qst_from_files <- function(trait_qst_file, neutral_dir, threshold_percentile,
                                   output_file, sanity_check = FALSE) {
  load(trait_qst_file)
  trait_result <- result
  trait_qst <- trait_result$qst

  trait_id_from_file <- sub("_trait_qst\\.RData$", "", basename(trait_qst_file))
  obs_stats_file <- sub("_trait_qst\\.RData$", "_obs_stats.RData", trait_qst_file)
  if (file.exists(obs_stats_file)) {
    load(obs_stats_file)
  } else {
    trait_meta <- list(trait_id = trait_id_from_file, chr = "unknown")
    obs_stats <- NULL
  }
  if (!exists("trait_meta")) {
    trait_meta <- list(trait_id = trait_id_from_file, chr = "unknown")
  }

  trait_both_neg <- FALSE
  if (!is.null(trait_result$both_var_negative) && !is.na(trait_result$both_var_negative) &&
      as.integer(trait_result$both_var_negative) == 1L) {
    trait_both_neg <- TRUE
  } else if (!is.null(obs_stats)) {
    trait_both_neg <- is_both_neg(obs_stats)
  }

  trait_abc_na <- FALSE
  if (!is.null(trait_result$abc_na) && as.integer(trait_result$abc_na) == 1L) {
    trait_abc_na <- TRUE
  } else if (!trait_both_neg && is.na(trait_qst)) {
    trait_abc_na <- TRUE
  }

  trait_num <- function(key) {
    if (is.null(trait_result[[key]])) return(NA_real_)
    v <- suppressWarnings(as.numeric(trait_result[[key]]))
    if (length(v) != 1L || is.na(v)) NA_real_ else v
  }

  single_files <- list.files(neutral_dir, pattern = "^neutral_[0-9]+\\.RData$", full.names = TRUE)
  batch_files <- list.files(neutral_dir, pattern = "^neutral_batch_[0-9]+\\.RData$", full.names = TRUE)

  neutral_qst_list <- list()
  neutral_fst_list <- list()
  neutral_both_neg_list <- list()
  neutral_abc_na_list <- list()
  idx_count <- 0

  for (f in single_files) {
    tryCatch({
      load(f)
      if (!is.null(result$fst)) {
        row_both_neg <- !is.null(result$both_var_negative) &&
          as.integer(result$both_var_negative) == 1L
        row_abc_na <- (!is.null(result$abc_na) && as.integer(result$abc_na) == 1L) ||
          (is.na(result$qst) && !row_both_neg)
        idx_count <- idx_count + 1
        neutral_qst_list[[idx_count]] <- result$qst
        neutral_fst_list[[idx_count]] <- result$fst
        neutral_both_neg_list[[idx_count]] <- row_both_neg
        neutral_abc_na_list[[idx_count]] <- row_abc_na
      }
    }, error = function(e) warning(paste("Failed to load single file:", f, "-", e$message)))
  }

  for (f in batch_files) {
    tryCatch({
      load(f)
      if (!is.null(result$results)) {
        df <- result$results
        row_both_neg <- if ("both_var_negative" %in% names(df)) {
          as.integer(df$both_var_negative) == 1L
        } else {
          rep(FALSE, nrow(df))
        }
        row_abc_na <- if ("abc_na" %in% names(df)) {
          as.integer(df$abc_na) == 1L
        } else {
          is.na(df$qst) & !row_both_neg
        }
        idx_count <- idx_count + 1
        neutral_qst_list[[idx_count]] <- df$qst
        neutral_fst_list[[idx_count]] <- df$fst
        neutral_both_neg_list[[idx_count]] <- row_both_neg
        neutral_abc_na_list[[idx_count]] <- row_abc_na
      }
    }, error = function(e) warning(paste("Failed to load batch file:", f, "-", e$message)))
  }

  if (length(neutral_qst_list) > 0) {
    neutral_qst_values <- unlist(neutral_qst_list, use.names = FALSE)
    neutral_fst_values <- unlist(neutral_fst_list, use.names = FALSE)
    neutral_both_neg <- unlist(neutral_both_neg_list, use.names = FALSE)
    neutral_abc_na <- unlist(neutral_abc_na_list, use.names = FALSE)
  } else {
    neutral_qst_values <- numeric(0)
    neutral_fst_values <- numeric(0)
    neutral_both_neg <- logical(0)
    neutral_abc_na <- logical(0)
  }

  valid_neutral <- !is.na(neutral_qst_values) & !neutral_both_neg & !neutral_abc_na
  n_valid_neutral <- sum(valid_neutral)
  if (n_valid_neutral < 2L) {
    threshold <- NA_real_
  } else {
    threshold <- quantile(neutral_qst_values[valid_neutral], threshold_percentile)
  }

  adaptive <- if (trait_both_neg || trait_abc_na || is.na(trait_qst) || is.na(threshold)) NA else trait_qst > threshold
  adaptive_str <- ifelse(is.na(adaptive), NA, ifelse(adaptive, "yes", "no"))

  if (sanity_check) {
    output <- data.frame(
      trait_id = trait_meta$trait_id,
      chr = trait_meta$chr,
      among_pop_sd = round(obs_num(obs_stats, "among_pop_sd"), 6),
      within_pop_sd = round(obs_num(obs_stats, "within_pop_sd"), 6),
      ext_sd = round(obs_num(obs_stats, "ext_sd"), 6),
      prior_QST = round(obs_num(obs_stats, "QST"), 6),
      QST = round(trait_qst, 6),
      QST_lo = round(trait_num("qst_lo"), 6),
      QST_hi = round(trait_num("qst_hi"), 6),
      abc_ratioVext = round(trait_num("abc_ratioVext"), 6),
      anova_ratioVext = round(if (!is.na(trait_num("anova_ratioVext"))) trait_num("anova_ratioVext") else obs_num(obs_stats, "ratioVext"), 6),
      ratioVext_source = if (!is.null(trait_result$ratioVext_source)) trait_result$ratioVext_source else "abc",
      both_var_negative = as.integer(trait_both_neg),
      abc_na = as.integer(trait_abc_na),
      threshold_percentile = threshold_percentile,
      threshold_value = round(as.numeric(threshold), 6),
      n_neutral = n_valid_neutral,
      adaptive = adaptive_str
    )
  } else {
    output <- data.frame(
      trait_id = trait_meta$trait_id,
      chr = trait_meta$chr,
      QST = round(trait_qst, 6),
      QST_lo = round(trait_num("qst_lo"), 6),
      QST_hi = round(trait_num("qst_hi"), 6),
      abc_ratioVext = round(trait_num("abc_ratioVext"), 6),
      anova_ratioVext = round(if (!is.na(trait_num("anova_ratioVext"))) trait_num("anova_ratioVext") else obs_num(obs_stats, "ratioVext"), 6),
      ratioVext_source = if (!is.null(trait_result$ratioVext_source)) trait_result$ratioVext_source else "abc",
      both_var_negative = as.integer(trait_both_neg),
      abc_na = as.integer(trait_abc_na),
      threshold_percentile = threshold_percentile,
      threshold_value = round(as.numeric(threshold), 6),
      adaptive = adaptive_str
    )
  }

  write.csv(output, output_file, row.names = FALSE, quote = FALSE)
  if (length(neutral_qst_values) > 0) {
    write.table(
      data.frame(fst = neutral_fst_values, qst = neutral_qst_values,
                 both_var_negative = as.integer(neutral_both_neg),
                 abc_na = as.integer(neutral_abc_na)),
      sub("\\.csv$", "_neutral_qst.txt", output_file),
      row.names = FALSE, col.names = TRUE, sep = "\t"
    )
  }
  invisible(output)
}

if (!interactive() && Sys.getenv("QST_ABC_SOURCE_ONLY", "") != "1") {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 4) {
    cat("Usage: Rscript aggregate_qst.R <trait_qst_file> <neutral_dir> <threshold_percentile> <output_file> [sanity_check]\n")
    quit(status = 1)
  }
  trait_qst_file <- args[1]
  neutral_dir <- args[2]
  threshold_percentile <- as.numeric(args[3])
  output_file <- args[4]
  sanity_check <- if (length(args) >= 5) as.logical(args[5]) else FALSE
  aggregate_qst_from_files(trait_qst_file, neutral_dir, threshold_percentile, output_file, sanity_check)
  cat("Results saved to:", output_file, "\n")
}
