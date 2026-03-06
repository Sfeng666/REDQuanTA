#!/usr/bin/env Rscript

# Plot Sample Structure Comparison Results
# Auto-discovers available sample structures from subdirectories.
#
# Usage:
#   Rscript plot_sample_structure_comparison.R <base_dir> <output_dir> <summary_stats>
#
# Arguments:
#   base_dir:       Directory containing n2_iX_rY/ subdirectories with TPR matrices
#   output_dir:     Where to save plots and tables
#   summary_stats:  Comma-separated summary statistics label (e.g. "QST,F_within_pop")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_sample_structure_comparison.R <base_dir> <output_dir> <summary_stats>")
}

base_dir <- args[1]
output_dir <- args[2]
summary_stats <- args[3]
summary_stats_str <- gsub(",", ".", summary_stats)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

theme_publication <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(size = base_size, color = "black"),
      axis.title = element_text(size = base_size + 1, face = "bold"),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.8),
      legend.key.size = unit(0.8, "lines")
    )
}

# Auto-discover sample structure subdirectories (n2_iX_rY pattern)
struct_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
struct_dirs <- struct_dirs[grepl("^n2_i\\d+_r\\d+$", struct_dirs)]

if (length(struct_dirs) == 0) {
  stop("No sample structure subdirectories (n2_iX_rY) found in: ", base_dir)
}

sp_structure <- data.frame(
  n_strain = as.integer(sub("n2_i(\\d+)_r\\d+", "\\1", struct_dirs)),
  n_replicate = as.integer(sub("n2_i\\d+_r(\\d+)", "\\1", struct_dirs)),
  stringsAsFactors = FALSE
)
sp_structure$struct_id <- struct_dirs
sp_structure$total_samples <- sp_structure$n_strain * sp_structure$n_replicate * 2
sp_structure <- sp_structure[order(sp_structure$total_samples, sp_structure$n_strain), ]

cat("Discovered", nrow(sp_structure), "sample structures:\n")
print(sp_structure[, c("struct_id", "n_strain", "n_replicate", "total_samples")])

load_data <- function(chr_name) {
  results_list <- list()
  for (i in 1:nrow(sp_structure)) {
    struct_id <- sp_structure$struct_id[i]
    result_file <- file.path(base_dir, struct_id, chr_name,
                             paste0("tpr_fpr_matrix_", chr_name, ".csv"))

    if (file.exists(result_file)) {
      cat("Reading:", result_file, "\n")
      data <- read.csv(result_file, header = TRUE, row.names = 1, check.names = FALSE)

      adaptive_rows <- grep("^QST_", rownames(data), value = TRUE)
      if (length(adaptive_rows) > 0) {
        adaptive_row <- adaptive_rows[length(adaptive_rows)]
        tpr_values <- as.numeric(data[adaptive_row, ])
        ratios <- colnames(data)

        for (j in seq_along(tpr_values)) {
          results_list[[length(results_list) + 1]] <- data.frame(
            n_strain = sp_structure$n_strain[i],
            n_replicate = sp_structure$n_replicate[i],
            total_samples = sp_structure$total_samples[i],
            ratio_venv_vtotal = ratios[j],
            TPR = tpr_values[j],
            chromosome = chr_name,
            stringsAsFactors = FALSE
          )
        }
      }
    } else {
      warning(paste("File not found:", result_file))
    }
  }
  return(results_list)
}

list_auto <- load_data("autosomes")
list_chrx <- load_data("chrX")
results_list <- c(list_auto, list_chrx)

if (length(results_list) == 0) {
  stop("No results found!")
}

df_results <- do.call(rbind, results_list)
df_results$ratio_venv_vtotal_num <- as.numeric(gsub("VEratio_", "", df_results$ratio_venv_vtotal))
df_results$n_replicate_factor <- factor(df_results$n_replicate,
                                         levels = sort(unique(df_results$n_replicate)))
df_results$chromosome <- factor(df_results$chromosome,
                                 levels = c("autosomes", "chrX"),
                                 labels = c("Autosomes", "X chromosome"))

cat("\nData summary:\n")
print(head(df_results))

df_results$ratio_label <- paste0("V[E]~'/'~V[G]~'='~'",
                                 format(df_results$ratio_venv_vtotal_num,
                                        scientific = TRUE, digits = 2), "'")

x_breaks <- sort(unique(df_results$total_samples))
n_replicates <- length(unique(df_results$n_replicate))
df_results$group_id <- interaction(df_results$n_replicate_factor, df_results$chromosome)

p_combined <- ggplot(df_results, aes(x = total_samples, y = TPR,
                                      color = n_replicate_factor,
                                      shape = n_replicate_factor,
                                      linetype = chromosome)) +
  geom_line(aes(group = group_id), linewidth = 1) +
  geom_point(size = 2.4, stroke = 0.42) +
  facet_wrap(~ratio_label, ncol = 3, labeller = label_parsed) +
  scale_x_continuous(breaks = x_breaks) +
  scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9,
                        name = "No. of\nreplicates") +
  scale_shape_manual(values = c(16, 17, 15, 18)[1:n_replicates],
                    name = "No. of\nreplicates") +
  scale_linetype_manual(values = c("Autosomes" = "solid", "X chromosome" = "dashed"),
                        name = "Chromosome") +
  guides(color = guide_legend(override.aes = list(linetype = 0, size = 3), order = 1),
         shape = guide_legend(override.aes = list(linetype = 0, size = 3), order = 1),
         linetype = guide_legend(order = 2, keywidth = unit(2.5, "lines"))) +
  labs(
    title = expression(bold("Power to detect adaptive "*Q[ST]*" across different sample structures")),
    subtitle = paste("Summary stats:", summary_stats),
    x = "Total sample size (No. of strains \u00d7 No. of replicates \u00d7 2 populations)",
    y = "True Positive Rate (Power)"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_publication(base_size = 10) +
  theme(
    legend.position = "right",
    panel.spacing = unit(0.8, "lines"),
    legend.box = "vertical",
    strip.text = element_text(size = 8)
  )

plot_file <- file.path(output_dir, paste0("power_comparison_combined_", summary_stats_str, ".pdf"))
ggsave(plot_file, p_combined, width = 170, height = 120, units = "mm", dpi = 300)
cat("Saved plot to:", plot_file, "\n")

df_results$samples_per_pop <- df_results$total_samples / 2
x_breaks_per_pop <- sort(unique(df_results$samples_per_pop))

p_per_pop <- ggplot(df_results, aes(x = samples_per_pop, y = TPR,
                                     color = n_replicate_factor,
                                     shape = n_replicate_factor,
                                     linetype = chromosome)) +
  geom_line(aes(group = group_id), linewidth = 1) +
  geom_point(size = 2.4, stroke = 0.42) +
  facet_wrap(~ratio_label, ncol = 3, labeller = label_parsed) +
  scale_x_continuous(breaks = x_breaks_per_pop) +
  scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9,
                        name = "No. of\nreplicates") +
  scale_shape_manual(values = c(16, 17, 15, 18)[1:n_replicates],
                    name = "No. of\nreplicates") +
  scale_linetype_manual(values = c("Autosomes" = "solid", "X chromosome" = "dashed"),
                        name = "Chromosome") +
  guides(color = guide_legend(override.aes = list(linetype = 0, size = 3), order = 1),
         shape = guide_legend(override.aes = list(linetype = 0, size = 3), order = 1),
         linetype = guide_legend(order = 2, keywidth = unit(2.5, "lines"))) +
  labs(
    title = expression(bold("Power to detect adaptive "*Q[ST]*" across different sample structures")),
    subtitle = paste("Summary stats:", summary_stats),
    x = "Sample size per population (No. of strains \u00d7 No. of replicates)",
    y = "True Positive Rate (Power)"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_publication(base_size = 10) +
  theme(
    legend.position = "right",
    panel.spacing = unit(0.8, "lines"),
    legend.box = "vertical",
    strip.text = element_text(size = 8)
  )

plot_file_per_pop <- file.path(output_dir, paste0("power_comparison_combined_per_pop_", summary_stats_str, ".pdf"))
ggsave(plot_file_per_pop, p_per_pop, width = 170, height = 120, units = "mm", dpi = 300)
cat("Saved per-pop plot to:", plot_file_per_pop, "\n")

for (chr_factor in unique(df_results$chromosome)) {
  chr_df <- df_results %>% filter(chromosome == chr_factor)

  chr_wide <- chr_df %>%
    select(n_strain, n_replicate, total_samples, ratio_venv_vtotal_num, TPR) %>%
    mutate(V_ratio = paste0("V_ratio_", as.character(ratio_venv_vtotal_num))) %>%
    select(-ratio_venv_vtotal_num) %>%
    pivot_wider(names_from = V_ratio, values_from = TPR) %>%
    arrange(total_samples, n_strain)

  chr_name <- if (chr_factor == "X chromosome") "chrX" else "autosomes"
  txt_file <- file.path(output_dir, paste0("power_summary_table_", chr_name, "_", summary_stats_str, ".txt"))
  write.table(chr_wide, file = txt_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Saved TPR table to:", txt_file, "\n")
}

cat("\nDone!\n")
