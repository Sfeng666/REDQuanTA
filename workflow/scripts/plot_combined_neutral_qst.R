#!/usr/bin/env Rscript
# Combine neutral QST/FST from trait batches (memory-safe) and plot pooled distributions.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  stop("Usage: Rscript plot_combined_neutral_qst.R <results_dir> [out_dir] [pdf_only]")
}
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
results_dir <- args[1]
out_dir <- if (length(args) >= 2) args[2] else file.path(results_dir, "plots")
pdf_only <- length(args) >= 3 && tolower(args[3]) %in% c("pdf_only", "pdf-only", "pdf")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(script_dir, "neutral_qst_memory_safe.R"))

qst_path <- file.path(results_dir, "summary", "qst_results_abc.csv")
if (!file.exists(qst_path)) {
  qst_path <- file.path(results_dir, "qst_results_abc.csv")
}
if (!file.exists(qst_path)) {
  stop("Trait QST table not found under ", results_dir)
}
validation_traits <- fread(qst_path, select = "trait_id")$trait_id

batch_files <- sort(list.files(
  results_dir,
  pattern = "^neutral_batch_[0-9]+\\.RData$",
  full.names = TRUE,
  recursive = TRUE
))
trait_from_path <- sub(".*/trait_([^/]+)/.*", "\\1", batch_files)
batch_files <- batch_files[trait_from_path %in% validation_traits]
if (length(batch_files) == 0) {
  stop("No neutral_batch_*.RData files for validation traits under ", results_dir)
}

cat("Streaming", length(batch_files), "neutral batch files (no full in-memory combine)...\n")
acc <- stream_neutral_batches(batch_files)
n_qst <- acc$n_pairs
n_traits <- length(validation_traits)

cat("Traits with neutral batches:", n_traits, "\n")
cat("Combined neutral (fst, qst) pairs:", format(n_qst, big.mark = ","), "\n")

save_neutral_accumulator(acc, out_dir)

summary_lines <- c(
  sprintf("n_traits_with_batches: %d", n_traits),
  sprintf("n_batch_files: %d", length(batch_files)),
  sprintf("n_fst_qst_pairs: %d", n_qst),
  sprintf("n_fst_lt0: %d", acc$n_fst_lt0),
  sprintf("expected_if_complete: %d (traits x 1000 neutral per trait)", n_traits * 1000L),
  sprintf("plot_method: full KDE (adjust=%s) + peak counter-shift + unit AUC; fst_bins=[%.2f,1]",
          Sys.getenv("NEUTRAL_PLOT_KDE_ADJUST", "1.05"), min(acc$fst_edges))
)
writeLines(summary_lines, file.path(out_dir, "combined_neutral_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

theme_pub <- function() {
  theme_classic(base_size = 11, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle = element_text(size = 10, colour = "grey30", hjust = 0),
      legend.position = c(0.72, 0.88),
      legend.background = element_rect(fill = alpha("white", 0.92), colour = NA),
      axis.line = element_line(linewidth = 0.5, colour = "black")
    )
}

x_label_fst_or_qst <- expression(italic(F)[ST] ~ "or" ~ italic(Q)[ST])

remove_deprecated_neutral_plot_outputs <- function(out_dir) {
  patterns <- c(
    "^neutral_qst_distribution_combined\\.",
    "^qst_distribution_publication(_pooled|_by_class)?\\.",
    "^qst_distribution_adaptive_only\\.",
    "^peak_decomposition_by_trait_group\\."
  )
  for (pat in patterns) {
    for (f in list.files(out_dir, pattern = pat, full.names = TRUE)) {
      unlink(f)
    }
  }
}

med_qst <- weighted_quantile_from_hist(acc$qst_counts, acc$edges, 0.5)
p95_qst <- weighted_quantile_from_hist(acc$qst_counts, acc$edges, 0.95)

trait_qst <- fread(qst_path)
trait_qst <- trait_qst[is.finite(QST), .(value = QST)]
n_trait <- nrow(trait_qst)
trait_counts <- hist_counts_in_edges(trait_qst$value, acc$edges)

# Full KDE + peak counter-shift on actual x, then unit AUC for comparability.
fst_fit <- hist_to_density_smooth(acc$fst_counts, acc$fst_edges)
qst_fit <- hist_to_density_smooth(acc$qst_counts, acc$edges)
trait_fit <- hist_to_density_smooth(trait_counts, acc$edges)

fst_curve <- normalize_density_auc(fst_fit$curve)
qst_curve <- normalize_density_auc(qst_fit$curve)
trait_curve <- normalize_density_auc(trait_fit$curve)

plot_dt <- rbindlist(list(
  fst_curve[, .(value = x, density = y, series = "Neutral FST (same draws)")],
  qst_curve[, .(value = x, density = y, series = "Neutral QST (ABC)")],
  trait_curve[, .(value = x, density = y, series = "Trait QST (ABC)")]
))
plot_dt[, series := factor(
  series,
  levels = c("Neutral FST (same draws)", "Neutral QST (ABC)", "Trait QST (ABC)")
)]

peak_meta <- data.table(
  series = c("Neutral FST (same draws)", "Neutral QST (ABC)", "Trait QST (ABC)"),
  hist_mode_x = c(fst_fit$mode_x, qst_fit$mode_x, trait_fit$mode_x),
  kde_peak_x = c(fst_fit$kde_peak_x, qst_fit$kde_peak_x, trait_fit$kde_peak_x),
  peak_shift = c(fst_fit$peak_shift, qst_fit$peak_shift, trait_fit$peak_shift),
  curve_peak_x = c(fst_fit$peak_x, qst_fit$peak_x, trait_fit$peak_x),
  kde_adjust = c(fst_fit$kde_adjust, qst_fit$kde_adjust, trait_fit$kde_adjust),
  n_obs = c(n_qst, n_qst, n_trait)
)
fwrite(peak_meta, file.path(out_dir, "neutral_fst_qst_trait_qst_overlay_peaks.csv"))

cols <- c(
  "Neutral FST (same draws)" = "#4E79A7",
  "Neutral QST (ABC)" = "#59A14F",
  "Trait QST (ABC)" = "#E15759"
)

trapz_area <- trapz_area_xy
fwrite(plot_dt, file.path(out_dir, "neutral_fst_qst_trait_qst_overlay_curves.csv"))
area_overlay <- plot_dt[, .(trapezoid_area = trapz_area(value, density)), by = series]
fwrite(area_overlay, file.path(out_dir, "neutral_fst_qst_trait_qst_overlay_areas.csv"))

p_overlay <- ggplot(plot_dt, aes(x = value, y = density, colour = series)) +
  geom_line(linewidth = 0.95) +
  scale_colour_manual(values = cols, name = NULL) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Neutral FST, neutral QST, and trait QST distributions",
    subtitle = sprintf(
      "Full KDE (adjust=1.05) + peak counter-shift, unit AUC; neutral n = %s (%s FST<0); trait n = %d",
      format(n_qst, big.mark = ","), format(acc$n_fst_lt0, big.mark = ","), n_trait
    ),
    x = x_label_fst_or_qst,
    y = "Density"
  ) +
  theme_pub()

ggsave(file.path(out_dir, "neutral_fst_qst_trait_qst_overlay.pdf"), p_overlay, width = 6.5, height = 4.2)
if (!pdf_only) {
  ggsave(file.path(out_dir, "neutral_fst_qst_trait_qst_overlay.png"), p_overlay, width = 6.5, height = 4.2, dpi = 600)
}

x_show <- 0.08
p_overlay_full <- ggplot(plot_dt, aes(x = value, y = density, colour = series)) +
  geom_line(linewidth = 0.95) +
  geom_vline(xintercept = c(0, 1), linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  scale_colour_manual(values = cols, name = NULL) +
  coord_cartesian(xlim = c(-x_show, 1 + x_show)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Neutral FST, neutral QST, and trait QST (full KDE, unit AUC)",
    subtitle = sprintf(
      "Same curves as overlay, wider x; KDE adjust=1.05, peak counter-shift; FST bins [%.2f,1] (n=%s FST<0); trait n = %d",
      min(acc$fst_edges), format(acc$n_fst_lt0, big.mark = ","), n_trait
    ),
    x = x_label_fst_or_qst,
    y = "Density"
  ) +
  theme_pub()

fwrite(plot_dt, file.path(out_dir, "neutral_fst_qst_trait_qst_overlay_full_kde_curves.csv"))
fwrite(area_overlay, file.path(out_dir, "neutral_fst_qst_trait_qst_overlay_full_kde_areas.csv"))
ggsave(file.path(out_dir, "neutral_fst_qst_trait_qst_overlay_full_kde.pdf"), p_overlay_full, width = 6.8, height = 4.3)
if (!pdf_only) {
  ggsave(file.path(out_dir, "neutral_fst_qst_trait_qst_overlay_full_kde.png"), p_overlay_full, width = 6.8, height = 4.3, dpi = 600)
}

# Drop superseded outputs if present from older runs.
for (f in list.files(out_dir, pattern = "^neutral_fst_qst_trait_qst_overlay_(full_kde_norm01|unit01)", full.names = TRUE)) {
  unlink(f)
}
remove_deprecated_neutral_plot_outputs(out_dir)

cat("Plots written to", out_dir, "\n")
