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
adaptive_qst <- as.numeric(strsplit(args[3], ",")[[1]])
ve_ratios <- as.numeric(strsplit(args[4], ",")[[1]])
threshold_percentile <- as.numeric(args[5])
output_file <- args[6]

cat("Aggregating performance evaluation results\n")
cat("Results directory:", results_dir, "\n")
cat("Chromosome type:", chr_type, "\n")
cat("Adaptive QST levels:", paste(adaptive_qst, collapse = ", "), "\n")
cat("V_E/V_G ratios:", paste(ve_ratios, collapse = ", "), "\n")
cat("Threshold percentile:", threshold_percentile, "\n")

#' Load QST estimates from RData files in a directory
#' @param dir_path Directory containing batch RData files
#' @param pattern Pattern to match files
#' @param combo Filter by specific summary stat combination (NULL = return all)
#' @return If combo is NULL and multi-combo data: data.frame with combo and qst columns
#'         If combo is specified or single-combo data: vector of QST estimates
load_qst_estimates <- function(dir_path, pattern = "*.RData", combo = NULL) {
  files <- list.files(dir_path, pattern = "\\.RData$", full.names = TRUE)
  if (length(files) == 0) {
    warning(paste("No RData files found in:", dir_path))
    return(numeric(0))
  }
  
  all_results <- list()
  is_multi_combo <- FALSE
  
  for (f in files) {
    tryCatch({
      load(f)  # Loads 'result'
      
      # Check if this is multi-combo mode
      if (!is.null(result$is_multi_combo) && result$is_multi_combo) {
        is_multi_combo <- TRUE
      }
      
      if (result$mode == "batch_neutral") {
        # Neutral batch
        if (is_multi_combo && "combo" %in% names(result$results)) {
          all_results[[length(all_results) + 1]] <- result$results[, c("combo", "qst")]
        } else {
          all_results[[length(all_results) + 1]] <- data.frame(
            combo = "default", qst = result$results$qst
          )
        }
      } else if (result$mode == "batch_evaluate") {
        # Adaptive batch
        if (is_multi_combo && "combo" %in% names(result$results)) {
          all_results[[length(all_results) + 1]] <- data.frame(
            combo = result$results$combo,
            qst = result$results$estimated_qst
          )
        } else {
          all_results[[length(all_results) + 1]] <- data.frame(
            combo = "default", qst = result$results$estimated_qst
          )
        }
      } else if (result$mode == "neutral") {
        all_results[[length(all_results) + 1]] <- data.frame(combo = "default", qst = result$qst)
      } else if (result$mode == "evaluate") {
        all_results[[length(all_results) + 1]] <- data.frame(combo = "default", qst = result$estimated_qst)
      }
    }, error = function(e) {
      warning(paste("Error loading file:", f, "-", e$message))
    })
  }
  
  if (length(all_results) == 0) {
    return(if (is_multi_combo && is.null(combo)) data.frame(combo = character(), qst = numeric()) else numeric(0))
  }
  
  combined <- do.call(rbind, all_results)
  combined <- combined[!is.na(combined$qst), ]
  
  if (!is.null(combo)) {
    # Filter to specific combination
    return(combined$qst[combined$combo == combo])
  } else if (is_multi_combo) {
    # Return all data with combo column
    return(combined)
  } else {
    # Return just the QST values
    return(combined$qst)
  }
}

#' Get unique combinations from multi-combo results
get_unique_combos <- function(dir_path) {
  files <- list.files(dir_path, pattern = "\\.RData$", full.names = TRUE)
  combos <- c()
  
  for (f in files) {
    tryCatch({
      load(f)
      if (!is.null(result$is_multi_combo) && result$is_multi_combo && 
          "combo" %in% names(result$results)) {
        combos <- unique(c(combos, result$results$combo))
      }
    }, error = function(e) {})
  }
  
  return(combos)
}

# Initialize results matrices
n_qst <- length(adaptive_qst)
n_ratios <- length(ve_ratios)

cat("\n=== Loading and processing results ===\n")

# Check if this is multi-combo mode by examining first neutral directory
first_neutral_dir <- file.path(results_dir, "neutral_ratio_1")
combo_list <- get_unique_combos(first_neutral_dir)
is_multi_combo <- length(combo_list) > 0

if (is_multi_combo) {
  cat("\nMulti-combo mode detected:", length(combo_list), "combinations\n")
} else {
  combo_list <- c("default")
  cat("\nSingle-combo mode\n")
}

# Process each combination
for (combo in combo_list) {
  combo_name <- if (combo == "default") "" else paste0("_", gsub(",", "_", combo))
  cat("\n=== Processing combination:", combo, "===\n")
  
  # TPR matrix: rows = adaptive QST levels, columns = V_E/V_G ratios
  tpr_matrix <- matrix(NA, nrow = n_qst, ncol = n_ratios)
  rownames(tpr_matrix) <- paste0("QST_", adaptive_qst)
  colnames(tpr_matrix) <- paste0("VEratio_", ve_ratios)
  
  # FPR matrix: one row (neutral), columns = V_E/V_G ratios
  fpr_matrix <- matrix(NA, nrow = 1, ncol = n_ratios)
  rownames(fpr_matrix) <- "neutral"
  colnames(fpr_matrix) <- paste0("VEratio_", ve_ratios)
  
  # Threshold matrix: neutral QST threshold at each V_E ratio
  threshold_matrix <- matrix(NA, nrow = 1, ncol = n_ratios)
  rownames(threshold_matrix) <- "threshold"
  colnames(threshold_matrix) <- paste0("VEratio_", ve_ratios)
  
  # Step 1: Load neutral QST estimates and calculate thresholds for each V_E ratio
  cat("\nProcessing neutral QST estimates...\n")
  
  for (ratio_idx in 1:n_ratios) {
    ve_ratio <- ve_ratios[ratio_idx]
    neutral_dir <- file.path(results_dir, paste0("neutral_ratio_", ratio_idx))
    
    if (!dir.exists(neutral_dir)) {
      cat("  Warning: Neutral directory not found:", neutral_dir, "\n")
      next
    }
    
    if (is_multi_combo) {
      neutral_qst <- load_qst_estimates(neutral_dir, combo = combo)
    } else {
      neutral_qst <- load_qst_estimates(neutral_dir)
    }
    n_neutral <- length(neutral_qst)
    
    if (n_neutral > 0) {
      # Calculate threshold at given percentile
      threshold <- quantile(neutral_qst, probs = threshold_percentile, na.rm = TRUE)
      threshold_matrix[1, ratio_idx] <- threshold
      
      # Calculate FPR (false positive rate)
      fpr <- sum(neutral_qst > threshold) / n_neutral
      fpr_matrix[1, ratio_idx] <- fpr
      
      cat("  V_E ratio", ratio_idx, "(", ve_ratio, "): n =", n_neutral, 
          ", threshold =", round(threshold, 4), ", FPR =", round(fpr, 4), "\n")
    } else {
      cat("  V_E ratio", ratio_idx, "(", ve_ratio, "): No data\n")
    }
  }
  
  # Step 2: Load adaptive QST estimates and calculate TPR for each condition
  cat("\nProcessing adaptive QST estimates...\n")
  
  for (qst_idx in 1:n_qst) {
    qst_value <- adaptive_qst[qst_idx]
    qst_str <- gsub("\\.", "_", sprintf("%.2f", qst_value))
    
    for (ratio_idx in 1:n_ratios) {
      ve_ratio <- ve_ratios[ratio_idx]
      adaptive_dir <- file.path(results_dir, paste0("adaptive_q", qst_str, "_r", ratio_idx))
      
      if (!dir.exists(adaptive_dir)) {
        cat("  Warning: Adaptive directory not found:", adaptive_dir, "\n")
        next
      }
      
      if (is_multi_combo) {
        adaptive_qst_estimates <- load_qst_estimates(adaptive_dir, combo = combo)
      } else {
        adaptive_qst_estimates <- load_qst_estimates(adaptive_dir)
      }
      n_adaptive <- length(adaptive_qst_estimates)
      
      if (n_adaptive > 0 && !is.na(threshold_matrix[1, ratio_idx])) {
        threshold <- threshold_matrix[1, ratio_idx]
        
        # Calculate TPR (true positive rate)
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
  
  # Output file for this combo
  if (is_multi_combo && combo != "default") {
    combo_output_file <- sub("\\.csv$", paste0(combo_name, ".csv"), output_file)
  } else {
    combo_output_file <- output_file
  }
  
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
        main = paste("True Positive Rate (TPR) -", chr_type, 
                     if(is_multi_combo) paste0("\nModel: ", combo) else ""))
  
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
  
  # Collect model ranking data
  if (!exists("model_ranking")) {
    model_ranking <- data.frame()
  }
  
  model_ranking <- rbind(model_ranking, data.frame(
    model = combo,
    mean_tpr = mean(tpr_matrix, na.rm = TRUE),
    mean_fpr = mean(fpr_matrix, na.rm = TRUE),
    n_stats = length(strsplit(combo, ",")[[1]])
  ))

}  # End of combo loop

# Step 5: Generate model ranking if multi-combo
if (is_multi_combo && nrow(model_ranking) > 1) {
  cat("\n=== Model Ranking ===\n")
  
  # Sort by mean TPR (descending)
  model_ranking <- model_ranking[order(-model_ranking$mean_tpr), ]
  model_ranking$rank <- 1:nrow(model_ranking)
  
  ranking_file <- sub("\\.csv$", "_model_ranking.csv", output_file)
  write.csv(model_ranking, ranking_file, row.names = FALSE)
  cat("Model ranking saved to:", ranking_file, "\n")
  
  # Print top 10
  cat("\nTop 10 models by mean TPR:\n")
  for (i in 1:min(10, nrow(model_ranking))) {
    row <- model_ranking[i, ]
    cat(sprintf("  %2d. %s (TPR=%.3f, FPR=%.3f, %d stats)\n",
                i, row$model, row$mean_tpr, row$mean_fpr, row$n_stats))
  }
} else {
  # Print summary for single combo
  cat("\n=== Summary ===\n")
  cat("Threshold percentile:", threshold_percentile, "\n")
  cat("Mean FPR across V_E ratios:", round(mean(fpr_matrix, na.rm = TRUE), 4), "\n")
  cat("Mean TPR across all conditions:", round(mean(tpr_matrix, na.rm = TRUE), 4), "\n")
  cat("\nTPR by adaptive QST level (averaged across V_E ratios):\n")
  for (i in 1:n_qst) {
    cat("  QST =", adaptive_qst[i], ":", round(mean(tpr_matrix[i, ], na.rm = TRUE), 4), "\n")
  }
}

cat("\nDone!\n")
