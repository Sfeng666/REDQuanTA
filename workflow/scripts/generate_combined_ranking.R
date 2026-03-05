#!/usr/bin/env Rscript
# Generate combined model ranking across all chromosomes
# Auto-discovers chromosome subdirectories under the results directory
#
# Usage:
#   Rscript generate_combined_ranking.R [results_dir]
#
# Arguments:
#   results_dir: Directory containing chromosome subdirectories (default: results/perf_eval)
#                Each subdirectory should contain tpr_fpr_matrix_*.csv files

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Resolve project root relative to this script's location
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  return(getwd())
}
script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, ".."))

# Default results directory (relative to project root)
default_results_dir <- file.path(project_root, "code", "chtc", "results", "perf_eval")
results_dir <- if (length(args) >= 1) normalizePath(args[1]) else default_results_dir

cat("Results directory:", results_dir, "\n")

# Auto-discover chromosome subdirectories
chr_dirs <- list.dirs(results_dir, recursive = FALSE, full.names = TRUE)
# Keep only directories that contain tpr_fpr CSV files (skip plots/, etc.)
chr_dirs <- chr_dirs[sapply(chr_dirs, function(d) {
  length(list.files(d, pattern = "^tpr_fpr_matrix_.*\\.csv$")) > 0
})]
chr_types <- basename(chr_dirs)

if (length(chr_types) == 0) {
  stop("No chromosome subdirectories with tpr_fpr_matrix CSV files found in: ", results_dir)
}
cat("Found chromosome types:", paste(chr_types, collapse = ", "), "\n")

# Function to extract TPR by V_E ratio from individual CSV files
extract_tpr_by_ratio <- function(chr_dir, chr_type) {
  # Get list of tpr_fpr CSV files (exclude model_ranking and heatmap files)
  pattern <- paste0("tpr_fpr_matrix_", chr_type, "_.*\\.csv$")
  files <- list.files(chr_dir, pattern = pattern, full.names = TRUE)
  files <- files[!grepl("model_ranking|heatmap", files)]

  all_results <- list()

  for (f in files) {
    tryCatch({
      df <- read.csv(f, stringsAsFactors = FALSE)

      # Extract model name from filename
      model_name <- sub(paste0("tpr_fpr_matrix_", chr_type, "_"), "", basename(f))
      model_name <- sub("\\.csv$", "", model_name)

      # Get TPR rows (skip threshold and FPR rows)
      tpr_rows <- df[grepl("^QST_", df$type), ]

      if (nrow(tpr_rows) > 0) {
        ratio_cols <- grep("^VEratio_", names(df), value = TRUE)
        mean_tpr <- sapply(ratio_cols, function(col) mean(tpr_rows[[col]], na.rm = TRUE))

        result <- data.frame(model = model_name, chr = chr_type, stringsAsFactors = FALSE)
        for (i in seq_along(ratio_cols)) {
          result[[paste0("TPR_r", i)]] <- mean_tpr[i]
        }
        all_results[[length(all_results) + 1]] <- result
      }
    }, error = function(e) {
      warning(paste("Error processing:", f, "-", e$message))
    })
  }

  if (length(all_results) > 0) do.call(rbind, all_results) else NULL
}

# Extract TPR data from each chromosome
all_chr_data <- list()
for (i in seq_along(chr_types)) {
  cat("Extracting TPR data from", chr_types[i], "...\n")
  df_chr <- extract_tpr_by_ratio(chr_dirs[i], chr_types[i])
  if (!is.null(df_chr)) {
    cat("  Found", nrow(df_chr), "models\n")
    all_chr_data[[length(all_chr_data) + 1]] <- df_chr
  } else {
    cat("  No models found\n")
  }
}

if (length(all_chr_data) == 0) stop("No TPR data found in any chromosome directory")

df_combined <- do.call(rbind, all_chr_data)

# Discover V_E ratio columns and their values
tpr_cols <- grep("^TPR_r", names(df_combined), value = TRUE)

# Get V_E ratio values from the column names in the source CSVs
# Read one file to discover the actual ratio values
sample_file <- list.files(chr_dirs[1], pattern = "^tpr_fpr_matrix_.*\\.csv$", full.names = TRUE)
sample_file <- sample_file[!grepl("model_ranking|heatmap", sample_file)][1]
sample_df <- read.csv(sample_file, stringsAsFactors = FALSE)
ratio_cols <- grep("^VEratio_", names(sample_df), value = TRUE)
ve_ratios <- as.numeric(sub("^VEratio_", "", ratio_cols))

cat("Calculating mean TPR across", length(chr_types), "chromosome types...\n")

# Aggregate by model (mean across chromosomes)
df_avg <- aggregate(
  df_combined[, tpr_cols, drop = FALSE],
  by = list(model = df_combined$model),
  FUN = mean, na.rm = TRUE
)

# Overall mean TPR across all ratios
df_avg$mean_TPR <- rowMeans(df_avg[, tpr_cols, drop = FALSE], na.rm = TRUE)

# Sort and rank
df_avg <- df_avg[order(-df_avg$mean_TPR), ]
df_avg$rank <- 1:nrow(df_avg)

# Rename TPR columns to descriptive names
for (i in seq_along(ve_ratios)) {
  old_name <- paste0("TPR_r", i)
  new_name <- paste0("TPR_VE_", ve_ratios[i])
  colnames(df_avg)[colnames(df_avg) == old_name] <- new_name
}

# Reorder columns
tpr_ve_cols <- paste0("TPR_VE_", ve_ratios)
df_avg <- df_avg[, c("rank", "model", tpr_ve_cols, "mean_TPR")]

# Round
for (col in c(tpr_ve_cols, "mean_TPR")) df_avg[[col]] <- round(df_avg[[col]], 4)

# Save
output_file <- if (length(args) >= 2) args[2] else file.path(results_dir, "combined_model_ranking.csv")
write.csv(df_avg, output_file, row.names = FALSE)
cat("\nCombined model ranking saved to:", output_file, "\n")

cat("\nTop 10 models by mean TPR across", paste(chr_types, collapse = " and "), ":\n")
cat(strrep("=", 60), "\n")
print(head(df_avg, 10))
cat("\nDone!\n")
