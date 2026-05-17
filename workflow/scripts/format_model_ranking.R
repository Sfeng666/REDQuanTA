#!/usr/bin/env Rscript
# Format combined model ranking for publication
# 1. Optional: write combined_model_ranking_publication.csv (same as input but publication stat names)
# 2. Keep only TPR at specified V_E/V_G ratio for the two-column table
# 3. Remove rows containing skewness or kurtosis
# 4. Rename summary stats to publication-friendly names
# 5. Output tab-delimited table and table legend
#
# Usage:
#   Rscript format_model_ranking.R [results_dir] [ve_ratio]
#
# Arguments:
#   results_dir: Directory containing combined_model_ranking.csv (default: code/chtc/results/perf_eval)
#   ve_ratio:    V_E/V_G ratio to select TPR column for (default: 1.0)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  return(getwd())
}
script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."))

default_results_dir <- file.path(project_root, "results", "perf_eval")
results_dir <- if (length(args) >= 1) normalizePath(args[1]) else default_results_dir
ve_ratio <- if (length(args) >= 2) as.numeric(args[2]) else 1.0

input_file <- file.path(results_dir, "combined_model_ranking.csv")
cat("Input:", input_file, "\n")
cat("V_E/V_G ratio:", ve_ratio, "\n")

df_full <- read.csv(input_file, stringsAsFactors = FALSE)

# Longer tokens first so e.g. ratioVextVtotal before ratioVext.
stat_map <- c(
  "among_pop_sd" = "SD_GB",
  "within_pop_sd" = "SD_GW",
  "ratioVbetweenVtotal" = "V_GB/(V_G + V_E)",
  "ratioVwithinVtotal" = "V_GW/(V_G + V_E)",
  "ratioVbetweenVext" = "V_GB/V_E",
  "ratioVextVtotal" = "V_E/(V_G + V_E)",
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

# Skewness / kurtosis filter
df_full <- df_full[!grepl("skewness_data|kurtosis_data", df_full$model), ]

# Publication-named full ranking (same columns and row order as filtered input)
df_pub <- df_full
df_pub$model <- vapply(df_pub$model, rename_model, character(1))
pub_file <- file.path(results_dir, "combined_model_ranking_publication.csv")
write.csv(df_pub, pub_file, row.names = FALSE)
cat("Publication combined ranking ->", pub_file, "\n")

tpr_col <- paste0("TPR_VE_", ve_ratio)
if (!tpr_col %in% names(df_pub)) stop("Column not found: ", tpr_col)

df <- df_pub[, c("model", tpr_col)]
colnames(df)[2] <- "TPR"

df <- df[order(-df$TPR), ]
df$TPR <- round(df$TPR, 4)

df <- df[, c("model", "TPR")]
colnames(df) <- c("Summary statistics", "True Positive Rate (TPR)")

output_file <- file.path(results_dir, "Table_model_ranking.txt")
write.table(df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Table saved to:", output_file, "\n")

cat("\nFormatted table:\n")
print(df, row.names = FALSE)

legend_text <- 'Table X. Ranking of ABC models by statistical power to detect adaptive trait differentiation.

Models are ranked by true positive rate (TPR) at V_E/V_G = 1.0 (i.e., extrinsic variance equals total additive genetic variance), averaged across 11 adaptive Q_ST levels (0.50-1.00) and both autosomes and X chromosome. Each model uses a distinct combination of summary statistics for ABC-based Q_ST estimation. False positive rate (FPR) is controlled at 0.05 across all models via dynamic thresholding. Summary statistic abbreviations: Q_ST, ratio of between-population genetic variance to total genetic variance [V_GB / (V_GB + 2 * V_GW)]; V_G = V_GB + V_GW (total additive genetic variance); V_E/V_G, ratio of extrinsic variance to V_G; SD_GB, standard deviation of between-population genetic variance; SD_GW, standard deviation of within-population genetic variance; SD_E, standard deviation of extrinsic variance; V_GB/V_GW, ratio of between-population to within-population genetic variance; V_GW/V_E, ratio of within-population genetic variance to extrinsic variance; V_GB/V_E, ratio of between-population genetic variance to extrinsic variance; V_GB/(V_G + V_E), V_GW/(V_G + V_E), and V_E/(V_G + V_E), ratios of V_GB, V_GW, or V_E to total phenotypic variance (V_G + V_E).'

legend_file <- file.path(results_dir, "Table_model_ranking_legend.txt")
writeLines(legend_text, legend_file)
cat("\nLegend saved to:", legend_file, "\n")
