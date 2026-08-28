#!/usr/bin/env Rscript
# ABC QST density for one trait-type collection.
# Usage: Rscript plot_qst_distribution.R <results_dir> [out_dir] [pdf_only] [title]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else {
  file.path(getwd(), "htcondor", "results", "traits_harmonizr_validation")
}
out_dir <- if (length(args) >= 2) args[2] else file.path(results_dir, "plots")
pdf_only <- length(args) >= 3 && tolower(args[3]) %in% c("pdf_only", "pdf-only", "pdf")
plot_title <- if (length(args) >= 4) args[4] else "Distribution of ABC-estimated QST"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

qst_path <- file.path(results_dir, "summary", "qst_results_abc.csv")
if (!file.exists(qst_path)) {
  qst_path <- file.path(results_dir, "qst_results_abc.csv")
}
if (!file.exists(qst_path)) {
  stop("Run aggregate_detection_results.py first: ", qst_path)
}
dt <- fread(qst_path)
dt[, adaptive := factor(adaptive, levels = c("no", "yes"))]
n <- nrow(dt)
n_ad <- sum(dt$adaptive == "yes", na.rm = TRUE)
pct_ad <- round(100 * n_ad / n, 1)

# --- Main density figure ---
p <- ggplot(dt, aes(x = QST, fill = adaptive, color = adaptive)) +
  geom_density(alpha = 0.35, linewidth = 0.9, adjust = 1.05) +
  geom_rug(aes(color = adaptive), alpha = 0.25, length = unit(0.02, "npc")) +
  scale_fill_manual(
    values = c("no" = "#4E79A7", "yes" = "#E15759"),
    labels = c("no" = "Not adaptive", "yes" = "Adaptive"),
    name = NULL
  ) +
  scale_color_manual(
    values = c("no" = "#2F5597", "yes" = "#B40426"),
    labels = c("no" = "Not adaptive", "yes" = "Adaptive"),
    name = NULL
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = plot_title,
    subtitle = sprintf("n = %d traits, %d adaptive (%.1f%%)", n, n_ad, pct_ad),
    x = expression(italic(Q)[ST] ~ "(ABC estimate)"),
    y = "Density"
  ) +
  theme_classic(base_size = 11, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 10, colour = "grey30", hjust = 0),
    legend.position = c(0.82, 0.88),
    legend.background = element_rect(fill = alpha("white", 0.9), colour = NA),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    axis.line = element_line(linewidth = 0.5, colour = "black")
  )

# --- Pooled: adaptive + non-adaptive on one curve ---
p_pooled <- ggplot(dt, aes(x = QST)) +
  geom_density(fill = "#4E79A7", colour = "#2F5597", alpha = 0.45, linewidth = 0.9, adjust = 1.05) +
  geom_rug(colour = "#2F5597", alpha = 0.2, length = unit(0.02, "npc")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Full distribution of ABC-estimated QST (all traits)",
    subtitle = sprintf(
      "Adaptive and non-adaptive combined; n = %d (%d adaptive, %.1f%%)",
      n, n_ad, pct_ad
    ),
    x = expression(italic(Q)[ST] ~ "(ABC estimate)"),
    y = "Density"
  ) +
  theme_classic(base_size = 11, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 10, colour = "grey30", hjust = 0),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    axis.line = element_line(linewidth = 0.5, colour = "black")
  )

ggsave(file.path(out_dir, "qst_distribution_publication_pooled.pdf"), p_pooled, width = 6.5, height = 4.2)
ggsave(file.path(out_dir, "qst_distribution_publication_by_class.pdf"), p, width = 6.5, height = 4.2)

# Main publication PDF: page 1 pooled, page 2 by adaptive class
pdf(file.path(out_dir, "qst_distribution_publication.pdf"), width = 6.5, height = 4.2)
print(p_pooled)
print(p)
dev.off()

if (!pdf_only) {
  ggsave(file.path(out_dir, "qst_distribution_publication_pooled.png"), p_pooled, width = 6.5, height = 4.2, dpi = 600)
  ggsave(file.path(out_dir, "qst_distribution_publication_by_class.png"), p, width = 6.5, height = 4.2, dpi = 600)
  png(file.path(out_dir, "qst_distribution_publication.png"), width = 6.5, height = 4.2, units = "in", res = 600)
  print(p_pooled)
  dev.off()
}

# --- Inset: adaptive only (if enough) ---
if (n_ad >= 5) {
  pad <- dt[adaptive == "yes"]
  p2 <- ggplot(pad, aes(x = QST)) +
    geom_histogram(bins = 20, fill = "#E15759", colour = "white", linewidth = 0.2) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title = sprintf("Adaptive traits only (n = %d)", n_ad),
      x = expression(italic(Q)[ST]),
      y = "Count"
    ) +
    theme_classic(base_size = 10)
  ggsave(file.path(out_dir, "qst_distribution_adaptive_only.pdf"), p2, width = 4.5, height = 3.2)
  if (!pdf_only) {
    ggsave(file.path(out_dir, "qst_distribution_adaptive_only.png"), p2, width = 4.5, height = 3.2, dpi = 300)
  }
}

cat("Publication QST plots written to", out_dir, "\n")
