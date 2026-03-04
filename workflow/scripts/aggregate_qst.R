#!/usr/bin/env Rscript
# Aggregate neutral QST results and compare to trait QST
#
# Usage:
#   Rscript aggregate_qst.R <trait_qst_file> <neutral_dir> <threshold_percentile> <output_file> [sanity_check]
#
# Example:
#   Rscript aggregate_qst.R results/trait_L0MQ04/L0MQ04_trait_qst.RData results/trait_L0MQ04 0.95 results/trait_L0MQ04/L0MQ04_result.csv FALSE
#
# Arguments:
#   1. trait_qst_file: Path to trait QST result RData
#   2. neutral_dir: Directory containing neutral QST RData files
#   3. threshold_percentile: Percentile for threshold (default: 0.95)
#   4. output_file: Path to save final result CSV
#   5. sanity_check: Whether to output detailed diagnostics (default: FALSE)

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

# Load trait QST result
load(trait_qst_file)
trait_result <- result
trait_qst <- trait_result$qst

# Also load trait metadata
obs_stats_file <- sub("_trait_qst\\.RData$", "_obs_stats.RData", trait_qst_file)
if (file.exists(obs_stats_file)) {
  load(obs_stats_file)
} else {
  trait_meta <- list(trait_id = "unknown", chr = "unknown")
  obs_stats <- trait_result$obs_stats
}

cat("Trait:", trait_meta$trait_id, "\n")
cat("Trait QST:", trait_qst, "\n")

# Collect all neutral QST values
# Supports both single-FST files (neutral_*.RData) and batch files (neutral_batch_*.RData)
single_files <- list.files(neutral_dir, pattern = "^neutral_[0-9]+\\.RData$", full.names = TRUE)
batch_files <- list.files(neutral_dir, pattern = "^neutral_batch_[0-9]+\\.RData$", full.names = TRUE)

cat("Found", length(single_files), "single neutral files\n")
cat("Found", length(batch_files), "batch neutral files\n")

neutral_qst_values <- c()
neutral_fst_values <- c()

# Process single-FST files
for (f in single_files) {
  tryCatch({
    load(f)
    if (!is.null(result$qst) && !is.na(result$qst)) {
      neutral_qst_values <- c(neutral_qst_values, result$qst)
      neutral_fst_values <- c(neutral_fst_values, result$fst)
    }
  }, error = function(e) {
    warning(paste("Failed to load single file:", f, "-", e$message))
  })
}

# Process batch files
for (f in batch_files) {
  tryCatch({
    load(f)
    if (!is.null(result$results)) {
      # Batch mode: results is a data frame with fst and qst columns
      valid <- !is.na(result$results$qst)
      neutral_qst_values <- c(neutral_qst_values, result$results$qst[valid])
      neutral_fst_values <- c(neutral_fst_values, result$results$fst[valid])
    }
  }, error = function(e) {
    warning(paste("Failed to load batch file:", f, "-", e$message))
  })
}

cat("Valid neutral QST values:", length(neutral_qst_values), "\n")

# Calculate threshold
threshold <- quantile(neutral_qst_values, threshold_percentile, na.rm = TRUE)
cat("Threshold (", threshold_percentile * 100, "%):", threshold, "\n")

# Determine if adaptive
adaptive <- if (is.na(trait_qst)) NA else trait_qst > threshold
adaptive_str <- ifelse(is.na(adaptive), NA, ifelse(adaptive, "yes", "no"))

cat("Adaptive:", adaptive_str, "\n")

# Create output
if (sanity_check) {
  # Detailed output with variance components
  output <- data.frame(
    trait_id = trait_meta$trait_id,
    chr = trait_meta$chr,
    among_pop_sd = round(obs_stats['among_pop_sd'], 6),
    within_pop_sd = round(obs_stats['within_pop_sd'], 6),
    ext_sd = round(obs_stats['ext_sd'], 6),
    prior_QST = round(obs_stats['QST'], 6),
    QST = round(trait_qst, 6),
    threshold_percentile = threshold_percentile,
    threshold_value = round(as.numeric(threshold), 6),
    n_neutral = length(neutral_qst_values),
    adaptive = adaptive_str
  )
  
  # Also save the full neutral distribution
  neutral_dist_file <- sub("\\.csv$", "_neutral_qst.txt", output_file)
  write.table(data.frame(fst = neutral_fst_values, qst = neutral_qst_values), 
              neutral_dist_file, row.names = FALSE, col.names = TRUE, sep = "\t")
  cat("Neutral distribution saved to:", neutral_dist_file, "\n")
  
} else {
  # Standard output
  output <- data.frame(
    trait_id = trait_meta$trait_id,
    chr = trait_meta$chr,
    QST = round(trait_qst, 6),
    threshold_percentile = threshold_percentile,
    threshold_value = round(as.numeric(threshold), 6),
    adaptive = adaptive_str
  )
}

write.csv(output, output_file, row.names = FALSE, quote = FALSE)
cat("Results saved to:", output_file, "\n")

