#!/usr/bin/env bash
# Optional post-processing: combined model ranking + publication table.
#
# Requires per-combo tpr_fpr_matrix_*_combo_*.csv under autosomes/ and chrX/
# (from run_aggregate_perf_eval_multicombo_fast.sh or equivalent aggregation).
#
# Writes under <perf_eval_root>:
#   combined_model_ranking.csv
#   combined_model_ranking_publication.csv
#   Table_model_ranking.txt
#   Table_model_ranking_legend.txt
#
# Usage:
#   bash run_perf_eval_publication_ranking.sh <perf_eval_root> [ve_ratio]
#
#   ve_ratio — V_E/V_G for the two-column Table_model_ranking.txt (default: 1.0)
#
# Environment:
#   RSCRIPT — path to Rscript (default: Rscript on PATH)
#   SKIP_COMBINED=1 — skip generate_combined_ranking.R (use existing combined CSV)
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT_BIN="${RSCRIPT:-$(command -v Rscript)}"
GENERATE_R="$SCRIPT_DIR/generate_combined_ranking.R"
FORMAT_R="$SCRIPT_DIR/format_model_ranking.R"

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <perf_eval_root> [ve_ratio]" >&2
  exit 1
fi

OUT_ROOT="$(cd "$1" && pwd)"
VE_RATIO="${2:-1.0}"

if [[ "${SKIP_COMBINED:-0}" != "1" ]]; then
  echo "Building combined_model_ranking.csv..."
  "$RSCRIPT_BIN" "$GENERATE_R" "$OUT_ROOT"
else
  echo "SKIP_COMBINED=1: using existing combined_model_ranking.csv"
fi

if [[ ! -f "$OUT_ROOT/combined_model_ranking.csv" ]]; then
  echo "ERROR: missing $OUT_ROOT/combined_model_ranking.csv" >&2
  exit 1
fi

echo "Formatting publication ranking (V_E/V_G = $VE_RATIO)..."
"$RSCRIPT_BIN" "$FORMAT_R" "$OUT_ROOT" "$VE_RATIO"

echo "Done. Key outputs:"
echo "  $OUT_ROOT/combined_model_ranking_publication.csv"
echo "  $OUT_ROOT/Table_model_ranking.txt"
