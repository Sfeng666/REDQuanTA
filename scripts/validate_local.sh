#!/bin/bash
# REDQuanTEA local smoke: Detection (2 traits) and Design (tiny grid).
set -euo pipefail

echo "============================================"
echo "REDQuanTEA local validation"
echo "============================================"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
cd "$REPO_DIR"

if ! command -v conda >/dev/null 2>&1; then
  if [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    # shellcheck source=/dev/null
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
  elif [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    # shellcheck source=/dev/null
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
  else
    echo "ERROR: conda not found."
    exit 1
  fi
fi

ENV_NAME="redquantea"
if conda env list | grep -qE '(^|[[:space:]])redquantea([[:space:]]|$)'; then
  ENV_NAME="redquantea"
elif conda env list | grep -qE '(^|[[:space:]])redquanta([[:space:]]|$)'; then
  ENV_NAME="redquanta"
  echo "Using existing conda env 'redquanta' (create 'redquantea' from environment.yml when you can)."
else
  echo "Creating redquantea conda environment..."
  conda env create -f environment.yml
  ENV_NAME="redquantea"
fi

echo ""
echo "Step 1: R packages ($ENV_NAME)"
conda run -n "$ENV_NAME" Rscript -e "library(abc); library(ggplot2); cat('R packages OK\n')"

echo ""
echo "Step 2: Snakemake version"
conda run -n "$ENV_NAME" snakemake --version

echo ""
echo "Step 3: Detection Module smoke (2 traits)"
conda run --no-capture-output -n "$ENV_NAME" snakemake \
    --snakefile workflow/Snakefile \
    --configfile config/config_detect_smoke.yaml \
    --directory "$REPO_DIR" \
    --cores 4

if [ ! -f results/detect_smoke/qst_results.csv ]; then
  echo "ERROR: missing results/detect_smoke/qst_results.csv"
  exit 1
fi
echo "Detection output:"
head -n 5 results/detect_smoke/qst_results.csv

echo ""
echo "Step 4: Design Module smoke (tiny grid)"
conda run --no-capture-output -n "$ENV_NAME" snakemake evaluate_all \
    --snakefile workflow/Snakefile \
    --configfile config/config_evaluate_smoke.yaml \
    --directory "$REPO_DIR" \
    --cores 4

if [ ! -f results/evaluate_smoke/combined_model_ranking.csv ]; then
  echo "ERROR: missing results/evaluate_smoke/combined_model_ranking.csv"
  exit 1
fi
ls -1 results/evaluate_smoke/*/tpr_fpr_matrix_*.csv
echo "Design ranking:"
cat results/evaluate_smoke/combined_model_ranking.csv

echo ""
echo "============================================"
echo "Validation PASSED"
echo "============================================"
