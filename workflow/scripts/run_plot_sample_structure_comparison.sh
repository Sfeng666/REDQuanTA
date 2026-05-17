#!/usr/bin/env bash
# Optional post-processing: power plots and TPR tables for sample-structure comparison.
#
# Run on any directory that contains n2_i{ind}_r{rep}/ subdirs with completed perf-eval
# (each struct has autosomes/ and chrX/ with tpr_fpr_matrix_*.csv).
#
# Usage:
#   bash run_plot_sample_structure_comparison.sh <results_base_dir> [summary_stats] [plots_dir]
#
# Examples:
#   bash workflow/scripts/run_plot_sample_structure_comparison.sh \
#     htcondor/results/validation_sample_struct_QST_ratioVbetweenVtotal
#
#   bash workflow/scripts/run_plot_sample_structure_comparison.sh \
#     results/sample_struct_comparison_QST_ratioVbetweenVtotal \
#     "QST,ratioVbetweenVtotal" results/sample_struct_comparison_QST_ratioVbetweenVtotal/plots
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT_BIN="${RSCRIPT:-$(command -v Rscript)}"

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <results_base_dir> [summary_stats] [plots_dir]" >&2
  exit 1
fi

BASE_DIR="$(cd "$1" && pwd)"
SUMMARY_STATS="${2:-QST,ratioVbetweenVtotal}"
PLOT_DIR="${3:-$BASE_DIR/plots}"

mkdir -p "$PLOT_DIR"
echo "Plotting sample-structure comparison"
echo "  base:   $BASE_DIR"
echo "  plots:  $PLOT_DIR"
echo "  stats:  $SUMMARY_STATS"

"$RSCRIPT_BIN" "$SCRIPT_DIR/plot_sample_structure_comparison.R" \
  "$BASE_DIR" "$PLOT_DIR" "$SUMMARY_STATS"

echo "Done. See $PLOT_DIR"
