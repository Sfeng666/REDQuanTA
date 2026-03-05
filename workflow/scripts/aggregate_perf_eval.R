#!/usr/bin/env Rscript
# Aggregate performance evaluation results (Module 2)
# Calculates TPR/FPR matrices and generates heatmaps
#
# Usage:
#   Rscript aggregate_perf_eval.R <results_dir> <chr_type> <adaptive_qst> <ve_ratios> <threshold> <output_file> [summary_stats]
#
# Example:
#   Rscript aggregate_perf_eval.R results/perf_eval/autosomes autosomes "0.50,0.55,0.60" "0.01,0.1,1.0" 0.95 results/perf_eval/autosomes/tpr_fpr_matrix_autosomes.csv
#
# Arguments:
#   1. results_dir: Directory containing neutral and adaptive results
#   2. chr_type: 'autosomes' or 'chrX'
#   3. adaptive_qst: Comma-separated adaptive QST levels
#   4. ve_ratios: Comma-separated V_E/V_G ratios
#   5. threshold_percentile: Percentile for threshold (default: 0.95)
#   6. output_file: Path to save TPR/FPR matrix CSV
#   7. summary_stats: (optional) Summary stats model name for model comparison

suppressPackageStartupMessages({
  library(grDevices)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  cat("Usage: Rscript aggregate_perf_eval.R <results_dir> <chr_type> <adaptive_qst> <ve_ratios> <threshold_percentile> <output_file>\n")
  quit(status = 1)
}

results_dir <- args[1]
chr_type <- args[2]
adaptive_qst_str <- strsplit(args[3], ",")[[1]]
ve_ratios_str <- strsplit(args[4], ",")[[1]]
adaptive_qst <- as.numeric(adaptive_qst_str)
ve_ratios <- as.numeric(ve_ratios_str)
threshold_percentile <- as.numeric(args[5])
output_file <- args[6]

cat("Aggregating performance evaluation results\n")
cat("Results directory:", results_dir, "\n")
cat("Chromosome type:", chr_type, "\n")
cat("Adaptive QST levels:", paste(adaptive_qst, collapse = ", "), "\n")
cat("V_E/V_G ratios:", paste(ve_ratios, collapse = ", "), "\n")
cat("Threshold percentile:", threshold_percentile, "\n")

#' Load QST estimates from RData files matching a glob pattern
#' Files are found in results_dir by pattern (flat naming, not subdirectories)
load_qst_from_files <- function(files) {
  if (length(files) == 0) return(numeric(0))
  
  all_qst <- c()
  for (f in files) {
    tryCatch({
      load(f)
      if (result$mode == "batch_neutral") {
        all_qst <- c(all_qst, result$results$qst)
      } else if (result$mode == "batch_evaluate") {
        all_qst <- c(all_qst, result$results$estimated_qst)
      } else if (result$mode == "neutral") {
        all_qst <- c(all_qst, result$qst)
      } else if (result$mode == "evaluate") {
        all_qst <- c(all_qst, result$estimated_qst)
      }
    }, error = function(e) {
      warning(paste("Error loading file:", f, "-", e$message))
    })
  }
  return(all_qst[!is.na(all_qst)])
}

n_qst <- length(adaptive_qst)
n_ratios <- length(ve_ratios)

cat("\n=== Loading and processing results ===\n")

combo_list <- c("default")

for (combo in combo_list) {
  tpr_matrix <- matrix(NA, nrow = n_qst, ncol = n_ratios)
  rownames(tpr_matrix) <- paste0("QST_", adaptive_qst)
  colnames(tpr_matrix) <- paste0("VEratio_", ve_ratios)

  fpr_matrix <- matrix(NA, nrow = 1, ncol = n_ratios)
  rownames(fpr_matrix) <- "neutral"
  colnames(fpr_matrix) <- paste0("VEratio_", ve_ratios)

  threshold_matrix <- matrix(NA, nrow = 1, ncol = n_ratios)
  rownames(threshold_matrix) <- "threshold"
  colnames(threshold_matrix) <- paste0("VEratio_", ve_ratios)

  cat("\nProcessing neutral QST estimates...\n")

  all_rdata <- list.files(results_dir, pattern = "\\.RData$", full.names = TRUE, recursive = TRUE)

  for (ratio_idx in 1:n_ratios) {
    ve_ratio <- ve_ratios[ratio_idx]
    prefix <- paste0("neutral_r", ve_ratios_str[ratio_idx], "_b")
    files <- all_rdata[grepl(prefix, basename(all_rdata), fixed = TRUE)]

    neutral_qst <- load_qst_from_files(files)
    n_neutral <- length(neutral_qst)

    if (n_neutral > 0) {
      threshold <- quantile(neutral_qst, probs = threshold_percentile, na.rm = TRUE)
      threshold_matrix[1, ratio_idx] <- threshold
      fpr <- sum(neutral_qst > threshold) / n_neutral
      fpr_matrix[1, ratio_idx] <- fpr

      cat("  V_E ratio", ratio_idx, "(", ve_ratio, "): n =", n_neutral,
          ", threshold =", round(threshold, 4), ", FPR =", round(fpr, 4), "\n")
    } else {
      cat("  V_E ratio", ratio_idx, "(", ve_ratio, "): No data (prefix:", prefix, ")\n")
    }
  }

  cat("\nProcessing adaptive QST estimates...\n")

  for (qst_idx in 1:n_qst) {
    qst_value <- adaptive_qst[qst_idx]

    for (ratio_idx in 1:n_ratios) {
      ve_ratio <- ve_ratios[ratio_idx]
      prefix <- paste0("adaptive_q", adaptive_qst_str[qst_idx], "_r", ve_ratios_str[ratio_idx], "_b")
      files <- all_rdata[grepl(prefix, basename(all_rdata), fixed = TRUE)]

      adaptive_qst_estimates <- load_qst_from_files(files)
      n_adaptive <- length(adaptive_qst_estimates)

      if (n_adaptive > 0 && !is.na(threshold_matrix[1, ratio_idx])) {
        threshold <- threshold_matrix[1, ratio_idx]
        tpr <- sum(adaptive_qst_estimates > threshold) / n_adaptive
        tpr_matrix[qst_idx, ratio_idx] <- tpr
      }
    }

    cat("  QST =", qst_value, ": TPR across ratios =",
        paste(round(tpr_matrix[qst_idx, ], 3), collapse = ", "), "\n")
  }

  # Step 3: Save results to CSV for this combo
  cat("\n=== Saving results for", combo, "===\n")
  
  # Create output directory if needed
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  combo_output_file <- output_file
  
  # Combine all results into a single data frame
  results_df <- data.frame(
    type = c("threshold", "FPR", rownames(tpr_matrix)),
    rbind(threshold_matrix, fpr_matrix, tpr_matrix)
  )
  
  write.csv(results_df, combo_output_file, row.names = FALSE)
  cat("TPR/FPR matrix saved to:", combo_output_file, "\n")
  
  # Step 4: Generate heatmap
  heatmap_file <- sub("\\.csv$", "_heatmap.pdf", combo_output_file)
  
  pdf(heatmap_file, width = 10, height = 8)
  
  # Prepare data for heatmap (only TPR, not threshold/FPR rows)
  heatmap_data <- tpr_matrix
  heatmap_data[is.na(heatmap_data)] <- 0
  
  # Create color palette (white to blue)
  colors <- colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100)
  
  # Plot heatmap
  par(mar = c(8, 8, 4, 2))
  image(1:ncol(heatmap_data), 1:nrow(heatmap_data), t(heatmap_data),
        col = colors, axes = FALSE,
        xlab = "", ylab = "",
        main = paste("True Positive Rate (TPR) -", chr_type))
  
  # Add axes
  axis(1, at = 1:ncol(heatmap_data), labels = ve_ratios, las = 2, cex.axis = 0.8)
  mtext("V_E / V_G Ratio", side = 1, line = 5)
  axis(2, at = 1:nrow(heatmap_data), labels = adaptive_qst, las = 2, cex.axis = 0.8)
  mtext("Adaptive QST", side = 2, line = 5)
  
  # Add values to cells
  for (i in 1:ncol(heatmap_data)) {
    for (j in 1:nrow(heatmap_data)) {
      val <- heatmap_data[j, i]
      text(i, j, sprintf("%.2f", val), cex = 0.7,
           col = ifelse(val > 0.5, "white", "black"))
    }
  }
  
  # Add color bar
  legend("topright", legend = c("0.0", "0.5", "1.0"),
         fill = c("white", "lightblue", "darkblue"),
         title = "TPR", cex = 0.8)
  
  dev.off()
  cat("Heatmap saved to:", heatmap_file, "\n")
  
}  # End of combo loop

cat("\n=== Summary ===\n")
cat("Threshold percentile:", threshold_percentile, "\n")
cat("Mean FPR across V_E ratios:", round(mean(fpr_matrix, na.rm = TRUE), 4), "\n")
cat("Mean TPR across all conditions:", round(mean(tpr_matrix, na.rm = TRUE), 4), "\n")
cat("\nTPR by adaptive QST level (averaged across V_E ratios):\n")
for (i in 1:n_qst) {
  cat("  QST =", adaptive_qst[i], ":", round(mean(tpr_matrix[i, ], na.rm = TRUE), 4), "\n")
}

cat("\nDone!\n")
