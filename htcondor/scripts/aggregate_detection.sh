#!/usr/bin/env bash
# Aggregate Detection Module per-trait CSVs and (optionally) plot one collection.
#
# Usage:
#   bash htcondor/scripts/aggregate_detection.sh RESULTS_DIR [SUMMARY_DIR]
#
# RESULTS_DIR contains trait_*/{id}_result.csv.
# Writes summary/qst_results_abc.csv and adaptive_trait_summary.csv.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/../.." && pwd)"
RESULTS_DIR="${1:?usage: aggregate_detection.sh RESULTS_DIR [SUMMARY_DIR]}"
SUMMARY_DIR="${2:-$RESULTS_DIR/summary}"
PLOT="${PLOT:-1}"
RSCRIPT="${RSCRIPT:-$(command -v Rscript || true)}"

python3 "$SCRIPT_DIR/aggregate_detection_results.py" \
  --results-dir "$RESULTS_DIR" \
  --summary-dir "$SUMMARY_DIR"

if [[ "$PLOT" == "1" && -n "$RSCRIPT" ]]; then
  mkdir -p "$RESULTS_DIR/plots"
  "$RSCRIPT" "$REPO/workflow/scripts/plot_qst_distribution.R" \
    "$RESULTS_DIR" "$RESULTS_DIR/plots" pdf_only
fi
