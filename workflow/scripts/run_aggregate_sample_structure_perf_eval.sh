#!/usr/bin/env bash
# Optional: run aggregate_perf_eval.R for every n2_i*_r* structure under a results tree.
#
# Use after HTCondor perf-eval DAGs finish (before plotting). Re-writes
# tpr_fpr_matrix_{autosomes,chrX}.csv and heatmaps under each structure.
#
# Usage:
#   bash run_aggregate_sample_structure_perf_eval.sh <results_base_dir> \
#     [adaptive_qst] [ve_ratios] [threshold]
#
# Defaults match validation_sample_struct_QST_ratioVbetweenVtotal:
#   adaptive_qst=0.8   ve_ratios=0.1,1.0,10.0   threshold=0.95
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT_BIN="${RSCRIPT:-$(command -v Rscript)}"
AGG_R="$SCRIPT_DIR/aggregate_perf_eval.R"

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <results_base_dir> [adaptive_qst] [ve_ratios] [threshold]" >&2
  exit 1
fi

BASE_DIR="$(cd "$1" && pwd)"
ADAPTIVE_QST="${2:-0.8}"
VE_RATIOS="${3:-0.1,1.0,10.0}"
THRESH="${4:-0.95}"

mapfile -t structs < <(find "$BASE_DIR" -maxdepth 1 -type d -name 'n2_i*_r*' | sort)
if [[ ${#structs[@]} -eq 0 ]]; then
  echo "ERROR: no n2_i*_r* directories under $BASE_DIR" >&2
  exit 1
fi

echo "Aggregating ${#structs[@]} structures under $BASE_DIR"
echo "  adaptive_qst=$ADAPTIVE_QST  ve_ratios=$VE_RATIOS  threshold=$THRESH"

for struct_dir in "${structs[@]}"; do
  sid="$(basename "$struct_dir")"
  for chr in autosomes chrX; do
    chr_dir="$struct_dir/$chr"
    [[ -d "$chr_dir" ]] || continue
    out="$chr_dir/tpr_fpr_matrix_${chr}.csv"
    echo "[$sid/$chr]"
    "$RSCRIPT_BIN" "$AGG_R" \
      "$chr_dir" "$chr" "$ADAPTIVE_QST" "$VE_RATIOS" "$THRESH" "$out"
  done
done

echo "Done. Re-run plotting with run_plot_sample_structure_comparison.sh"
