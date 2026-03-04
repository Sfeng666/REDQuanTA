#!/usr/bin/env Rscript
# Format combined model ranking for publication
# 1. Keep only TPR at specified V_E/V_G ratio
# 2. Remove rows containing skewness or kurtosis
#    (ABC estimation based on them may assume specific distributions of
#     simulated trait data that don't necessarily agree with empirical data)
# 3. Rename summary stats to publication-friendly names
# 4. Output tab-delimited table and table legend
#
# Usage:
#   Rscript format_model_ranking.R [results_dir] [ve_ratio]
#
# Arguments:
#   results_dir: Directory containing combined_model_ranking.csv (default: results/perf_eval)
#   ve_ratio:    V_E/V_G ratio to select TPR column for (default: 1.0)

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

# Defaults
default_results_dir <- file.path(project_root, "code", "chtc", "results", "perf_eval")
results_dir <- if (length(args) >= 1) normalizePath(args[1]) else default_results_dir
ve_ratio <- if (length(args) >= 2) as.numeric(args[2]) else 1.0

input_file <- file.path(results_dir, "combined_model_ranking.csv")
cat("Input:", input_file, "\n")
cat("V_E/V_G ratio:", ve_ratio, "\n")

df <- read.csv(input_file, stringsAsFactors = FALSE)

# 1. Keep only TPR at specified V_E/V_G ratio
tpr_col <- paste0("TPR_VE_", ve_ratio)
if (!tpr_col %in% names(df)) stop("Column not found: ", tpr_col)
df <- df[, c("model", tpr_col)]
colnames(df)[2] <- "TPR"

# 2. Remove rows with skewness_data or kurtosis_data
df <- df[!grepl("skewness_data|kurtosis_data", df$model), ]

# 3. Rename summary stats to publication-friendly names
stat_map <- c(
  "among_pop_sd" = "SD_GB",
  "within_pop_sd" = "SD_GW",
  "ext_sd" = "SD_E",
  "ratioVext" = "V_E/V_G",
  "F_among_pop" = "V_GB/V_GW",
  "F_within_pop" = "V_GW/V_E",
  "QST" = "Q_ST"
)

stat_names_sorted <- names(stat_map)[order(-nchar(names(stat_map)))]

rename_model <- function(model_str) {
  s <- model_str
  for (i in 1:5) {
    for (nm in stat_names_sorted) {
      s <- sub(nm, paste0("<<", stat_map[nm], ">>"), s, fixed = TRUE)
    }
  }
  parts <- regmatches(s, gregexpr("<<[^>]+>>", s))[[1]]
  parts <- gsub("<<|>>", "", parts)
  paste(parts, collapse = ", ")
}

df$model <- sapply(df$model, rename_model)

# Sort by TPR descending and round
df <- df[order(-df$TPR), ]
df$TPR <- round(df$TPR, 4)

# Final table: two columns only
df <- df[, c("model", "TPR")]
colnames(df) <- c("Summary statistics", "True Positive Rate (TPR)")

# Write tab-delimited table
output_file <- file.path(results_dir, "Table_model_ranking.txt")
write.table(df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Table saved to:", output_file, "\n")

cat("\nFormatted table:\n")
print(df, row.names = FALSE)

# 4. Write table legend
legend_text <- 'Table X. Ranking of ABC models by statistical power to detect adaptive trait differentiation.

Models are ranked by true positive rate (TPR) at V_E/V_G = 1.0 (i.e., extrinsic variance equals total additive genetic variance), averaged across 11 adaptive Q_ST levels (0.50-1.00) and both autosomes and X chromosome. Each model uses a distinct combination of summary statistics for ABC-based Q_ST estimation. False positive rate (FPR) is controlled at 0.05 across all models via dynamic thresholding. Summary statistic abbreviations: Q_ST, ratio of between-population genetic variance to total genetic variance [V_GB / (V_GB + 2 * V_GW)]; V_E/V_G, ratio of extrinsic variance to total genetic variance; SD_GB, standard deviation of between-population genetic variance; SD_GW, standard deviation of within-population genetic variance; SD_E, standard deviation of extrinsic variance; V_GB/V_GW, ratio of between-population to within-population genetic variance; V_GW/V_E, ratio of within-population genetic variance to extrinsic variance.'

legend_file <- file.path(results_dir, "Table_model_ranking_legend.txt")
writeLines(legend_text, legend_file)
cat("\nLegend saved to:", legend_file, "\n")
