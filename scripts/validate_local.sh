#!/bin/bash
# REDQuanTA Local Validation Script
# Run this to validate the workflow without Docker

set -e

echo "============================================"
echo "REDQuanTA Local Validation"
echo "============================================"

# Check conda environment
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found. Please install Miniforge or Miniconda."
    exit 1
fi

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"

cd "$REPO_DIR"

echo ""
echo "Step 1: Checking environment..."
echo "-------------------------------------------"

# Check if redquanta environment exists
if ! conda env list | grep -q "redquanta"; then
    echo "Creating redquanta conda environment..."
    conda env create -f environment.yml
fi

# Activate environment (this won't persist but we'll use conda run)
echo "Environment 'redquanta' found."

echo ""
echo "Step 2: Verifying R packages..."
echo "-------------------------------------------"
conda run -n redquanta Rscript -e "library(abc); library(ggplot2); cat('R packages OK\n')"

echo ""
echo "Step 3: Verifying Snakemake..."
echo "-------------------------------------------"
conda run -n redquanta snakemake --version

echo ""
echo "Step 4: Running Module 1 dry-run..."
echo "-------------------------------------------"
conda run -n redquanta snakemake \
    --snakefile workflow/Snakefile \
    --configfile config/config_detect.yaml \
    --directory "$REPO_DIR" \
    --cores 1 \
    -n

echo ""
echo "Step 5: Running Module 2 dry-run..."
echo "-------------------------------------------"
conda run -n redquanta snakemake evaluate_all \
    --snakefile workflow/Snakefile \
    --configfile config/config_evaluate.yaml \
    --directory "$REPO_DIR" \
    --cores 1 \
    -n

echo ""
echo "============================================"
echo "Validation PASSED"
echo "============================================"
echo ""
echo "To run the full workflow:"
echo "  conda activate redquanta"
echo "  snakemake --configfile config/config_detect.yaml --cores 4"
echo ""
