#!/usr/bin/env Rscript
# Plot TPR performance for performance evaluation (Module 2)
# Generates line plots showing TPR across adaptive QST levels and V_E/V_G ratios
# Uses cowplot for publication-ready two-panel layout
#
# Arguments:
#   1. results_dir: Directory containing TPR/FPR matrix CSVs
#   2. output_dir: Directory for plot outputs
#   3. summary_stats_models: (optional) Semicolon-separated list of summary stat model names
#                           If multiple, will generate model comparison/ranking

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(reshape2)
  library(cowplot)
})

# Publication-ready theme for Wiley/MER journals
# Requirements: min 9pt font, 300 dpi, color-blind friendly
theme_publication <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(size = base_size, color = "black"),
      axis.title = element_text(size = base_size + 1, face = "bold"),
      plot.title = element_text(size = base_size, hjust = 0.5),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0.5),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(0.3, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.8),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(1.0, "lines"),
      plot.margin = margin(5, 3, 3, 3)
    )
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript plot_tpr_performance.R <results_dir> <output_dir> [summary_stats_models]\n")
  cat("\nExamples:\n")
  cat("  # Single model:\n")
  cat("  Rscript plot_tpr_performance.R results/perf_eval plots/\n")
  cat("\n  # Multiple models (for ranking):\n")
  cat("  Rscript plot_tpr_performance.R results/perf_eval plots/ 'QST,ratioVext;QST;ratioVext'\n")
  quit(status = 1)
}

results_dir <- args[1]
output_dir <- args[2]
summary_stats_models <- if (length(args) >= 3) strsplit(args[3], ";")[[1]] else NULL

cat("Plotting TPR performance\n")
cat("Results directory:", results_dir, "\n")
cat("Output directory:", output_dir, "\n")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#' Load TPR/FPR matrix from CSV
#' @param file_path Path to CSV file
#' @return Data frame with TPR values
load_tpr_matrix <- function(file_path) {
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  return(df)
}

#' Prepare data for plotting
#' @param df Data frame from TPR/FPR matrix
#' @return Reshaped data frame for ggplot
prepare_plot_data <- function(df) {
  # Extract TPR rows (those starting with "QST_")
  tpr_rows <- df[grepl("^QST_", df$type), ]
  
  # Extract adaptive QST levels from type column
  tpr_rows$adaptive_qst <- as.numeric(gsub("QST_", "", tpr_rows$type))
  
  # Get V_E ratio columns
  ve_cols <- grep("^VEratio_", colnames(df), value = TRUE)
  
  # Reshape to long format
  plot_data <- melt(tpr_rows, 
                    id.vars = c("type", "adaptive_qst"),
                    measure.vars = ve_cols,
                    variable.name = "ve_ratio",
                    value.name = "TPR")
  
  # Extract numeric V_E ratio
  plot_data$ve_ratio_num <- as.numeric(gsub("VEratio_", "", plot_data$ve_ratio))
  plot_data$ve_ratio <- factor(plot_data$ve_ratio_num)
  
  return(plot_data)
}

#' Create TPR line plot for a single chromosome
#' @param plot_data Prepared plot data
#' @param chr_type Chromosome type
#' @param output_path Output file path
#' @param show_legend Whether to show legend
#' @param show_y_axis Whether to show y-axis labels
create_tpr_plot <- function(plot_data, chr_type, output_path = NULL, 
                            show_legend = TRUE, show_y_axis = TRUE) {
  
  chr_label <- ifelse(chr_type == "chrX", "X chromosome", "Autosomes")
  
  p <- ggplot(plot_data, aes(x = adaptive_qst, y = TPR, 
                              color = ve_ratio, group = ve_ratio)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.5, stroke = 0.3) +
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9,
                          name = expression(bold(V[E] / V[G]))) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = chr_label,
      x = expression(bold("Adaptive " * Q[ST] * " level")),
      y = if (show_y_axis) "True Positive Rate (Power)" else NULL
    ) +
    theme_publication(base_size = 10)
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  if (!show_y_axis) {
    p <- p + theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  
  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 170, height = 140, units = "mm", dpi = 300)
    cat("Saved:", output_path, "\n")
  }
  
  return(p)
}

#' Create combined two-panel plot using cowplot
#' @param data_list List of plot data for each chromosome
#' @param output_path Output file path
create_combined_plot <- function(data_list, output_path) {
  
  # Create left panel (Autosomes) - without legend
  p_autosomes <- create_tpr_plot(data_list[["autosomes"]], "autosomes", 
                                  show_legend = FALSE, show_y_axis = TRUE)
  
  # Create right panel (X chromosome) - with legend
  p_chrX <- create_tpr_plot(data_list[["chrX"]], "chrX", 
                             show_legend = TRUE, show_y_axis = FALSE)
  
  # Extract legend from chrX plot
  legend <- get_legend(p_chrX)
  
  # Remove legend from chrX plot for combining
  p_chrX_no_legend <- p_chrX + theme(legend.position = "none")
  
  # Create main title
  title <- ggdraw() + 
    draw_label(
      expression(bold("Power to detect adaptive " * Q[ST])),
      fontface = 'bold', size = 11, x = 0.5, hjust = 0.5
    )
  
  # Combine the two plots (autosomes gets more width to compensate for Y-axis labels)
  plots_row <- plot_grid(p_autosomes, p_chrX_no_legend, 
                          nrow = 1, rel_widths = c(1.18, 1), 
                          align = "h", axis = "tb")
  
  # Add legend to the right
  plots_with_legend <- plot_grid(plots_row, legend, 
                                  nrow = 1, rel_widths = c(2.2, 0.5))
  
  # Add title on top
  p_combined <- plot_grid(title, plots_with_legend, 
                           ncol = 1, rel_heights = c(0.08, 1))
  
  # Save combined plot - 180mm width for 'Large' figure per Wiley guidelines
  ggsave(output_path, p_combined, width = 180, height = 80, units = "mm", dpi = 300)
  cat("Saved:", output_path, "\n")
  
  return(p_combined)
}

#' Calculate mean TPR for model ranking
#' @param df TPR/FPR matrix data frame
#' @param model_name Name of the summary stats model
#' @return Named list with model name and mean TPR
calc_mean_tpr <- function(df, model_name) {
  # Extract TPR rows
  tpr_rows <- df[grepl("^QST_", df$type), ]
  ve_cols <- grep("^VEratio_", colnames(df), value = TRUE)
  
  # Calculate mean TPR across all conditions
  all_tpr <- as.matrix(tpr_rows[, ve_cols])
  mean_tpr <- mean(all_tpr, na.rm = TRUE)
  
  return(list(model = model_name, mean_tpr = mean_tpr))
}

# Process each chromosome type
chr_types <- c("autosomes", "chrX")

all_plot_data <- list()

for (chr_type in chr_types) {
  tpr_file <- file.path(results_dir, chr_type, paste0("tpr_fpr_matrix_", chr_type, ".csv"))
  
  if (!file.exists(tpr_file)) {
    cat("Warning: TPR file not found:", tpr_file, "\n")
    next
  }
  
  cat("\nProcessing", chr_type, "...\n")
  
  # Load and prepare data
  df <- load_tpr_matrix(tpr_file)
  plot_data <- prepare_plot_data(df)
  plot_data$chromosome <- chr_type
  
  all_plot_data[[chr_type]] <- plot_data
  
  # Create individual chromosome plot
  output_path <- file.path(output_dir, paste0("TPR_plot_", chr_type, ".pdf"))
  create_tpr_plot(plot_data, chr_type, output_path)
}

# Create combined two-panel plot if we have data for both chromosomes
if (length(all_plot_data) == 2) {
  cat("\nCreating combined two-panel plot...\n")
  
  output_path <- file.path(output_dir, "TPR_plot_combined.pdf")
  create_combined_plot(all_plot_data, output_path)
}

# Model ranking (if multiple models specified)
if (!is.null(summary_stats_models) && length(summary_stats_models) > 1) {
  cat("\nGenerating model ranking...\n")
  
  ranking_results <- list()
  
  for (model in summary_stats_models) {
    model_clean <- gsub(",", ".", model)
    model_results <- list()
    
    for (chr_type in chr_types) {
      # Look for model-specific result file
      tpr_file <- file.path(results_dir, model_clean, chr_type, 
                            paste0("tpr_fpr_matrix_", chr_type, ".csv"))
      
      if (file.exists(tpr_file)) {
        df <- load_tpr_matrix(tpr_file)
        result <- calc_mean_tpr(df, model)
        model_results[[chr_type]] <- result$mean_tpr
      }
    }
    
    if (length(model_results) > 0) {
      ranking_results[[model]] <- mean(unlist(model_results))
    }
  }
  
  if (length(ranking_results) > 0) {
    # Create ranking data frame
    ranking_df <- data.frame(
      model = names(ranking_results),
      mean_tpr = unlist(ranking_results)
    )
    ranking_df <- ranking_df[order(ranking_df$mean_tpr, decreasing = TRUE), ]
    ranking_df$rank <- 1:nrow(ranking_df)
    
    # Save ranking
    ranking_file <- file.path(output_dir, "model_ranking_by_mean_TPR.csv")
    write.csv(ranking_df, ranking_file, row.names = FALSE)
    cat("Model ranking saved to:", ranking_file, "\n")
    
    # Print ranking
    cat("\n=== Model Ranking by Mean TPR ===\n")
    for (i in 1:nrow(ranking_df)) {
      cat(sprintf("%d. %s: %.4f\n", ranking_df$rank[i], ranking_df$model[i], ranking_df$mean_tpr[i]))
    }
  }
}

cat("\nDone!\n")
