#!/usr/bin/env bash
# Design Module: generate HTCondor DAGs for one summary-stat combo file and
# optionally submit them.
#
# Required:
#   COMBO_FILE    tab-separated combo file (one combo per line, e.g. QST<TAB>ratioVbetweenVtotal)
#   OUTPUT_DIR    results root (writes autosomes/ and chrX/ plus perf_eval_*.dag)
#
# Optional:
#   FST_INPUT_DIR     default: $REPO/data/example
#   NUM_REPEATS       default: 10000
#   NUM_NEUTRAL       default: 10000
#   NUM_SIM           default: 100000
#   BATCH_SIZE        default: 1000
#   CHR               autosomes|chrX|both (default: both)
#   ADAPTIVE_QST      comma-separated QST levels (default: DAG generator 0.50..1.00)
#   VE_RATIOS         comma-separated V_E/V_G (default: 0.01,0.1,1,10,100)
#   FLOOR_POLICY      default: ridge_floor
#   FLOOR_ALPHA       default: 0.1
#   JOB_PRIORITY      default: 1000
#   DRY_RUN           1 = generate DAGs only
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/../.." && pwd)"
PREP="$SCRIPT_DIR/prepare_design_dag.py"

COMBO_FILE="${COMBO_FILE:?set COMBO_FILE}"
OUTPUT_DIR="${OUTPUT_DIR:?set OUTPUT_DIR}"
FST_INPUT_DIR="${FST_INPUT_DIR:-$REPO/data/example}"
NUM_REPEATS="${NUM_REPEATS:-10000}"
NUM_NEUTRAL="${NUM_NEUTRAL:-10000}"
NUM_SIM="${NUM_SIM:-100000}"
BATCH_SIZE="${BATCH_SIZE:-1000}"
CHR="${CHR:-both}"
ADAPTIVE_QST="${ADAPTIVE_QST:-}"
VE_RATIOS="${VE_RATIOS:-}"
FLOOR_POLICY="${FLOOR_POLICY:-ridge_floor}"
FLOOR_ALPHA="${FLOOR_ALPHA:-0.1}"
JOB_PRIORITY="${JOB_PRIORITY:-1000}"
DRY_RUN="${DRY_RUN:-0}"
CODE_DIR="${CODE_DIR:-$REPO/workflow/scripts}"
R_ENV_TARBALL="${R_ENV_TARBALL:-}"
SUB_NAME="${SUB_NAME:-perf_eval_fast.sub}"
REQUEST_MEMORY="${REQUEST_MEMORY:-4GB}"

if [[ -z "$R_ENV_TARBALL" ]]; then
  R_ENV_TARBALL="$(cd "$SCRIPT_DIR" && python3 -c 'from r_env_paths import resolve_r_env_tarball; print(resolve_r_env_tarball())')"
fi

mkdir -p "$OUTPUT_DIR"

prep_args=(
  --summary-stats "$COMBO_FILE"
  --output-dir "$OUTPUT_DIR"
  --chr "$CHR"
  --num-repeats "$NUM_REPEATS"
  --num-neutral "$NUM_NEUTRAL"
  --num-sim "$NUM_SIM"
  --batch-size "$BATCH_SIZE"
  --code-dir "$CODE_DIR"
  --fst-input-dir "$FST_INPUT_DIR"
  --sub-name "$SUB_NAME"
  --request-memory "$REQUEST_MEMORY"
  --floor-policy "$FLOOR_POLICY"
  --floor-alpha "$FLOOR_ALPHA"
  --r-env-tarball "$R_ENV_TARBALL"
)
if [[ -n "$ADAPTIVE_QST" ]]; then
  prep_args+=(--adaptive-qst "$ADAPTIVE_QST")
fi
if [[ -n "$VE_RATIOS" ]]; then
  prep_args+=(--ve-ratios "$VE_RATIOS")
fi

echo "Generating Design Module DAGs..."
python3 "$PREP" "${prep_args[@]}"

if [[ "$DRY_RUN" == "1" ]]; then
  echo "DRY_RUN=1; DAGs in $OUTPUT_DIR"
  ls -1 "$OUTPUT_DIR"/perf_eval_*.dag
  exit 0
fi

chroms=()
if [[ "$CHR" == "both" ]]; then
  chroms=(autosomes chrX)
else
  chroms=("$CHR")
fi

clusters=()
for chr in "${chroms[@]}"; do
  dag="$OUTPUT_DIR/perf_eval_${chr}.dag"
  [[ -f "$dag" ]] || { echo "ERROR: missing $dag" >&2; exit 1; }
  echo "Submitting $(basename "$dag") at priority $JOB_PRIORITY"
  out=$(condor_submit_dag -Force -Priority "$JOB_PRIORITY" "$dag")
  echo "$out"
  cluster=$(echo "$out" | sed -n 's/.*cluster \([0-9][0-9]*\).*/\1/p' | head -1)
  clusters+=("${cluster:-}")
done
echo "Submitted clusters: ${clusters[*]}"
printf '%s\n' "${clusters[*]}" > "$OUTPUT_DIR/clusters.txt"
