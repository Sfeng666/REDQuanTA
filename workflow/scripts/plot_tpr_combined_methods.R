#!/usr/bin/env Rscript
# TPR plot: same style as reference TPR_plot_combined.pdf
# - ABC lines from perf_eval combo 'QST, F_within_pop'
# - ANOVA and ANOVA without replication from reference eval_performance_detect_adaptive_qst/data/output
# Outputs: combined (2-panel), autosomes-only, chrX-only; full ratios and subset (0.1, 1, 10).

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(reshape2)
  library(cowplot)
})

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

# V_E/V_G ratios: 5 values matching reference (1e-02 to 1e+02)
RATIOS_FULL <- c("1e-02", "1e-01", "1e+00", "1e+01", "1e+02")
RATIOS_SUBSET <- c("1e-01", "1e+00", "1e+01")
RATIOS_SUBSET_LABELS <- c("0.1", "1", "10")

# Paths: assume run from repo root or from code/
find_paths <- function() {
  wd <- getwd()
  if (dir.exists(file.path(wd, "htcondor"))) {
    code_dir <- wd
    repo_root <- dirname(wd)
  } else {
    repo_root <- wd
    code_dir <- file.path(repo_root, "code")
  }
  perf_eval <- file.path(code_dir, "chtc", "results", "perf_eval")
  ref_out <- file.path(code_dir, "reference", "eval_performance_detect_adaptive_qst", "data", "output")
  list(perf_eval = perf_eval, ref_out = ref_out)
}
paths <- find_paths()
PERF_EVAL_DIR <- paths$perf_eval
REF_OUTPUT <- paths$ref_out
ANOVA_DIR <- file.path(REF_OUTPUT, "ANOVA_estimate", "generate_simulated_params_sumstats_pairedpriorsampling10-3to10_uniformprior_sampleF1_samplestr3x6_simx100000")
ANOVA_NOREP_DIR <- file.path(REF_OUTPUT, "ANOVA_estimate_withoutrep", "generate_simulated_params_sumstats_pairedpriorsampling10-3to10_uniformprior_sampleF1_samplestr3x6_simx100000_withoutrep")

get_tpr_df_anova <- function(path_matrix, v_ratios) {
  df <- read.table(path_matrix, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  df <- df[!grepl("mean_adaptive|neutral", rownames(df)), ]
  df <- df[, colnames(df) %in% v_ratios, drop = FALSE]
  df$adaptive_QST_level <- as.numeric(gsub("adaptive_", "", rownames(df)))
  rownames(df) <- NULL
  df
}

get_tpr_df_abc <- function(path_csv, v_ratios) {
  # ABC CSV: type, VEratio_0.01, VEratio_0.1, VEratio_1, VEratio_10, VEratio_100
  # Map to 1e-02, 1e-01, 1e+00, 1e+01, 1e+02
  map_abc_to_sci <- c(
    "VEratio_0.01" = "1e-02", "VEratio_0.1" = "1e-01", "VEratio_1" = "1e+00",
    "VEratio_10" = "1e+01", "VEratio_100" = "1e+02"
  )
  df <- read.csv(path_csv, stringsAsFactors = FALSE)
  tpr_rows <- df[grepl("^QST_", df$type), ]
  adaptive <- as.numeric(gsub("QST_", "", tpr_rows$type))
  out <- data.frame(adaptive_QST_level = adaptive)
  for (col in names(map_abc_to_sci)) {
    sci <- map_abc_to_sci[col]
    if (col %in% colnames(tpr_rows) && sci %in% v_ratios)
      out[[sci]] <- tpr_rows[[col]]
  }
  out <- out[, c("adaptive_QST_level", intersect(v_ratios, names(out))), drop = FALSE]
  out
}

get_plot_data_one_chr <- function(chr, v_ratios) {
  chr_file <- if (chr == "autosomes") "autosomes" else "chrX"
  # ABC: QST, F_within_pop
  abc_path <- file.path(PERF_EVAL_DIR, chr_file, paste0("tpr_fpr_matrix_", chr_file, "_QST_F_within_pop.csv"))
  # ANOVA
  anova_path <- file.path(ANOVA_DIR, paste0("matrix_tpr_fpr_", chr_file, ".txt"))
  anova_norep_path <- file.path(ANOVA_NOREP_DIR, paste0("matrix_tpr_fpr_", chr_file, ".txt"))

  dfs <- list()
  if (file.exists(abc_path)) {
    d_abc <- get_tpr_df_abc(abc_path, v_ratios)
    d_abc$method <- "ABC"
    dfs <- c(dfs, list(d_abc))
  }
  if (file.exists(anova_path)) {
    d_anova <- get_tpr_df_anova(anova_path, v_ratios)
    d_anova$method <- "ANOVA"
    dfs <- c(dfs, list(d_anova))
  }
  if (file.exists(anova_norep_path)) {
    d_norep <- get_tpr_df_anova(anova_norep_path, v_ratios)
    d_norep$method <- "ANOVA\nwithout replication"
    dfs <- c(dfs, list(d_norep))
  }
  if (length(dfs) == 0) return(NULL)
  combined <- do.call(rbind, dfs)
  plot_cols <- setdiff(colnames(combined), c("adaptive_QST_level", "method"))
  long <- melt(combined,
    id.vars = c("adaptive_QST_level", "method"),
    measure.vars = plot_cols,
    variable.name = "V_env_ratio", value.name = "TPR"
  )
  long$method <- factor(long$method, levels = c("ABC", "ANOVA", "ANOVA\nwithout replication"))
  long$V_env_ratio <- factor(as.character(long$V_env_ratio), levels = v_ratios)
  long$chromosome <- chr
  long
}

create_combined_plot <- function(df_auto, df_x, v_ratios, ratio_labels, output_path, title_suffix = "") {
  # Match line_plot_TPR_performance.Rmd: same expressions, linewidth 0.7, point size 1.5, stroke 0.3
  main_title <- expression(bold("Power to detect adaptive " * Q[ST] * " across methods and levels of extrinsic variance"))
  # title_suffix unused for combined; main_title is expression for draw_label

  p_autosomes <- ggplot(df_auto, aes(x = adaptive_QST_level, y = TPR,
    color = V_env_ratio, linetype = method, group = interaction(method, V_env_ratio))) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.5, stroke = 0.3) +
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9, labels = ratio_labels,
      name = expression(V[E] * " / " * V[G])) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = expression(Q[ST] * " method")) +
    labs(title = "Autosomes", x = NULL, y = "True Positive Rate (Power)") +
    ylim(0, 1) +
    theme_publication(base_size = 10) +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      plot.margin = margin(5, 5, 2, 5))

  p_chrX <- ggplot(df_x, aes(x = adaptive_QST_level, y = TPR,
    color = V_env_ratio, linetype = method, group = interaction(method, V_env_ratio))) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.5, stroke = 0.3) +
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9, labels = ratio_labels,
      name = expression(V[E] * " / " * V[G])) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = expression(Q[ST] * " method")) +
    guides(linetype = guide_legend(keywidth = unit(1.5, "lines"))) +
    labs(title = "X chromosome", x = expression(bold("Adaptive " * Q[ST] * " level")),
      y = "True Positive Rate (Power)") +
    ylim(0, 1) +
    theme_publication(base_size = 10) +
    theme(legend.position = "right", legend.box = "vertical", plot.margin = margin(2, 5, 5, 5))

  legend <- get_legend(p_chrX)
  p_chrX_noleg <- p_chrX + theme(legend.position = "none")
  title <- ggdraw() + draw_label(main_title, fontface = "bold", size = 11, x = 0.5, hjust = 0.5)
  plots_col <- plot_grid(p_autosomes, p_chrX_noleg, ncol = 1, align = "v", axis = "lr")
  plots_with_legend <- plot_grid(plots_col, legend, nrow = 1, rel_widths = c(1, 0.28))
  p_combined <- plot_grid(title, plots_with_legend, ncol = 1, rel_heights = c(0.05, 1))
  ggsave(output_path, p_combined, width = 160, height = 180, units = "mm", dpi = 300)
  invisible(p_combined)
}

create_single_plot <- function(df, chr_label, v_ratios, ratio_labels, output_path) {
  # Match line_plot_TPR_performance.Rmd section 6: linewidth 0.8, point 1.8, stroke 0.4, base_size 11
  # Legend names bold; linetype legend "Q_ST estimation method"; title/subtitle/x as in Rmd
  p <- ggplot(df, aes(x = adaptive_QST_level, y = TPR, color = V_env_ratio, linetype = method,
    group = interaction(method, V_env_ratio))) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8, stroke = 0.4) +
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9, labels = ratio_labels,
      name = expression(bold(V[E] * " / " * V[G]))) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted"),
      name = expression(bold(Q[ST] * " estimation method"))) +
    labs(
      title = expression(bold("Power to detect adaptive " * Q[ST] * " across methods and levels of extrinsic variance")),
      subtitle = paste0("Chromosome: ", chr_label),
      x = expression(bold("Adaptive " * Q[ST] * " level")),
      y = "True Positive Rate (Power)"
    ) +
    ylim(0, 1) +
    theme_publication(base_size = 11)
  ggsave(output_path, p, width = 170, height = 140, units = "mm", dpi = 300)
  invisible(p)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  output_dir <- REF_OUTPUT
  if (length(args) >= 1) output_dir <- args[1]
  if (length(args) >= 2) PERF_EVAL_DIR <- args[2]
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Full ratios
  df_auto <- get_plot_data_one_chr("autosomes", RATIOS_FULL)
  df_x <- get_plot_data_one_chr("chrX", RATIOS_FULL)
  if (is.null(df_auto) || is.null(df_x)) {
    message("Missing ABC or ANOVA data; check PERF_EVAL_DIR and REF_OUTPUT paths.")
    quit(status = 1)
  }
  create_combined_plot(df_auto, df_x, RATIOS_FULL, RATIOS_FULL,
    file.path(output_dir, "TPR_plot_combined_ABC_QST_Fwithin_pop.pdf"))
  create_single_plot(df_auto, "Autosomes", RATIOS_FULL, RATIOS_FULL,
    file.path(output_dir, "TPR_plot_autosomes_ABC_QST_Fwithin_pop.pdf"))
  create_single_plot(df_x, "X chromosome", RATIOS_FULL, RATIOS_FULL,
    file.path(output_dir, "TPR_plot_chrX_ABC_QST_Fwithin_pop.pdf"))

  # Subset: 0.1, 1, 10 only
  df_auto_s <- get_plot_data_one_chr("autosomes", RATIOS_SUBSET)
  df_x_s <- get_plot_data_one_chr("chrX", RATIOS_SUBSET)
  create_combined_plot(df_auto_s, df_x_s, RATIOS_SUBSET, RATIOS_SUBSET_LABELS,
    file.path(output_dir, "TPR_plot_combined_ABC_QST_Fwithin_pop_VEVG_subset.pdf"))
  create_single_plot(df_auto_s, "Autosomes", RATIOS_SUBSET, RATIOS_SUBSET_LABELS,
    file.path(output_dir, "TPR_plot_autosomes_ABC_QST_Fwithin_pop_VEVG_subset.pdf"))
  create_single_plot(df_x_s, "X chromosome", RATIOS_SUBSET, RATIOS_SUBSET_LABELS,
    file.path(output_dir, "TPR_plot_chrX_ABC_QST_Fwithin_pop_VEVG_subset.pdf"))

  message("Plots saved to ", output_dir)
}

main()
