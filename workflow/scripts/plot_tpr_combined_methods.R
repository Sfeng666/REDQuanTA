#!/usr/bin/env Rscript
# TPR plot: ABC (Design Module REDQuanTEA) vs aligned ANOVA / ANOVA without replication.
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

RATIOS_FULL <- c("1e-02", "1e-01", "1e+00", "1e+01", "1e+02")
RATIOS_FULL_LABELS <- c("0.01", "0.1", "1", "10", "100")
RATIOS_SUBSET <- c("1e-01", "1e+00", "1e+01")
RATIOS_SUBSET_LABELS <- c("0.1", "1", "10")

MAP_ABC_TO_SCI <- c(
  "VEratio_0.01" = "1e-02", "VEratio_0.1" = "1e-01", "VEratio_1" = "1e+00",
  "VEratio_10" = "1e+01", "VEratio_100" = "1e+02"
)

get_tpr_df_csv <- function(path_csv, v_ratios) {
  df <- read.csv(path_csv, stringsAsFactors = FALSE)
  tpr_rows <- df[grepl("^QST_", df$type), ]
  adaptive <- as.numeric(gsub("QST_", "", tpr_rows$type))
  out <- data.frame(adaptive_QST_level = adaptive)
  for (col in names(MAP_ABC_TO_SCI)) {
    sci <- MAP_ABC_TO_SCI[col]
    if (col %in% colnames(tpr_rows) && sci %in% v_ratios) {
      out[[sci]] <- tpr_rows[[col]]
    }
  }
  out[, c("adaptive_QST_level", intersect(v_ratios, names(out))), drop = FALSE]
}

resolve_paths <- function(args) {
  wd <- getwd()
  plots_dir <- if (length(args) >= 1) args[1] else file.path(wd, "plots")
  abc_perf_dir <- if (length(args) >= 2) args[2] else file.path(wd, "module2_original")
  anova_perf_dir <- if (length(args) >= 3) args[3] else file.path(abc_perf_dir, "anova")
  anova_norep_perf_dir <- if (length(args) >= 4) args[4] else file.path(abc_perf_dir, "anova_norep")
  abc_combo_suffix <- if (length(args) >= 5) args[5] else "QST_F_within_pop"
  output_tag <- if (length(args) >= 6) args[6] else "QST_Fwithin_pop"
  detector_suffix <- trimws(Sys.getenv("PERF_EVAL_DETECTOR_SUFFIX", ""))
  list(
    plots_dir = plots_dir,
    abc_perf_dir = abc_perf_dir,
    anova_perf_dir = anova_perf_dir,
    anova_norep_perf_dir = anova_norep_perf_dir,
    abc_combo_suffix = abc_combo_suffix,
    output_tag = output_tag,
    detector_suffix = detector_suffix
  )
}

resolve_abc_matrix_path <- function(chr_file, paths) {
  det <- paths$detector_suffix
  suffix <- paths$abc_combo_suffix
  prefix <- if (nzchar(det)) {
    paste0("tpr_fpr_matrix_", chr_file, det, "_")
  } else {
    paste0("tpr_fpr_matrix_", chr_file, "_")
  }
  chr_dir <- file.path(paths$abc_perf_dir, chr_file)

  # Prefer ranking lookup: "best", model name, or already a combo_XXXX id
  ranking <- file.path(chr_dir, paste0(prefix, "model_ranking.csv"))
  if (file.exists(ranking)) {
    rk <- read.csv(ranking, stringsAsFactors = FALSE)
    combo_id <- NULL
    if (identical(suffix, "best") || !nzchar(suffix)) {
      combo_id <- rk$model_id[1]
    } else if (grepl("^combo_[0-9]+$", suffix)) {
      combo_id <- suffix
    } else {
      hit <- which(rk$model == suffix | rk$model == gsub("_", ",", suffix, fixed = TRUE))
      if (length(hit)) combo_id <- rk$model_id[hit[1]]
    }
    if (!is.null(combo_id)) {
      candidate <- file.path(chr_dir, paste0(prefix, combo_id, ".csv"))
      if (file.exists(candidate)) return(candidate)
    }
  }

  direct <- file.path(chr_dir, paste0(prefix, suffix, ".csv"))
  if (file.exists(direct)) return(direct)
  for (alt in c("QST_F_within_pop", "QST,F_within_pop", "combo_0001")) {
    candidate <- file.path(chr_dir, paste0(prefix, alt, ".csv"))
    if (file.exists(candidate)) return(candidate)
  }
  NA_character_
}

get_plot_data_one_chr <- function(chr, v_ratios, paths) {
  chr_file <- if (chr == "autosomes") "autosomes" else "chrX"
  abc_path <- resolve_abc_matrix_path(chr_file, paths)
  anova_path <- file.path(paths$anova_perf_dir, chr_file, paste0("tpr_fpr_matrix_", chr_file, "_ANOVA.csv"))
  anova_norep_path <- file.path(paths$anova_norep_perf_dir, chr_file, paste0("tpr_fpr_matrix_", chr_file, "_ANOVA_norep.csv"))

  dfs <- list()
  if (!is.na(abc_path) && file.exists(abc_path)) {
    d_abc <- get_tpr_df_csv(abc_path, v_ratios)
    # Always label as ABC so it matches the plot legend factor levels
    # (ABC vs ANOVA vs ANOVA without replication).
    d_abc$method <- "ABC"
    message("ABC matrix (", chr_file, "): ", basename(abc_path))
    dfs <- c(dfs, list(d_abc))
  } else {
    warning("ABC TPR matrix not found for ", chr_file)
  }
  if (file.exists(anova_path)) {
    d_anova <- get_tpr_df_csv(anova_path, v_ratios)
    d_anova$method <- "ANOVA"
    dfs <- c(dfs, list(d_anova))
  }
  if (file.exists(anova_norep_path)) {
    d_norep <- get_tpr_df_csv(anova_norep_path, v_ratios)
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

create_combined_plot <- function(df_auto, df_x, v_ratios, ratio_labels, output_path) {
  main_title <- expression(bold("Power to detect adaptive " * Q[ST] * " across methods and levels of extrinsic variance"))

  p_autosomes <- ggplot(df_auto, aes(x = adaptive_QST_level, y = TPR,
    color = V_env_ratio, linetype = method, group = interaction(method, V_env_ratio))) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.5, stroke = 0.3) +
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9, labels = ratio_labels,
      name = expression(V[E] * " / " * V[G])) +
    scale_linetype_manual(
      values = c("ABC" = "solid", "ANOVA" = "dashed", "ANOVA\nwithout replication" = "dotted"),
      name = expression(Q[ST] * " method")) +
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
    scale_linetype_manual(
      values = c("ABC" = "solid", "ANOVA" = "dashed", "ANOVA\nwithout replication" = "dotted"),
      name = expression(Q[ST] * " method")) +
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
  p <- ggplot(df, aes(x = adaptive_QST_level, y = TPR, color = V_env_ratio, linetype = method,
    group = interaction(method, V_env_ratio))) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8, stroke = 0.4) +
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9, labels = ratio_labels,
      name = expression(bold(V[E] * " / " * V[G]))) +
    scale_linetype_manual(
      values = c("ABC" = "solid", "ANOVA" = "dashed", "ANOVA\nwithout replication" = "dotted"),
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
  paths <- resolve_paths(args)
  dir.create(paths$plots_dir, showWarnings = FALSE, recursive = TRUE)

  df_auto <- get_plot_data_one_chr("autosomes", RATIOS_FULL, paths)
  df_x <- get_plot_data_one_chr("chrX", RATIOS_FULL, paths)
  if (is.null(df_auto) || is.null(df_x)) {
    message("Missing ABC or aligned ANOVA data.")
    message("ABC dir: ", paths$abc_perf_dir)
    message("ANOVA dir: ", paths$anova_perf_dir)
    message("ANOVA norep dir: ", paths$anova_norep_perf_dir)
    quit(status = 1)
  }

  tag <- paths$output_tag
  create_combined_plot(df_auto, df_x, RATIOS_FULL, RATIOS_FULL_LABELS,
    file.path(paths$plots_dir, paste0("TPR_plot_combined_ABC_", tag, ".pdf")))
  create_single_plot(df_auto, "Autosomes", RATIOS_FULL, RATIOS_FULL_LABELS,
    file.path(paths$plots_dir, paste0("TPR_plot_autosomes_ABC_", tag, ".pdf")))
  create_single_plot(df_x, "X chromosome", RATIOS_FULL, RATIOS_FULL_LABELS,
    file.path(paths$plots_dir, paste0("TPR_plot_chrX_ABC_", tag, ".pdf")))

  df_auto_s <- get_plot_data_one_chr("autosomes", RATIOS_SUBSET, paths)
  df_x_s <- get_plot_data_one_chr("chrX", RATIOS_SUBSET, paths)
  create_combined_plot(df_auto_s, df_x_s, RATIOS_SUBSET, RATIOS_SUBSET_LABELS,
    file.path(paths$plots_dir, paste0("TPR_plot_combined_ABC_", tag, "_VEVG_subset.pdf")))
  create_single_plot(df_auto_s, "Autosomes", RATIOS_SUBSET, RATIOS_SUBSET_LABELS,
    file.path(paths$plots_dir, paste0("TPR_plot_autosomes_ABC_", tag, "_VEVG_subset.pdf")))
  create_single_plot(df_x_s, "X chromosome", RATIOS_SUBSET, RATIOS_SUBSET_LABELS,
    file.path(paths$plots_dir, paste0("TPR_plot_chrX_ABC_", tag, "_VEVG_subset.pdf")))

  message("Plots saved to ", paths$plots_dir)
}

main()
