# REDQuanTA Validation Guide

This document explains how to validate the REDQuanTA workflow in different environments.

## Validation Status

**Note**: Docker is not available on the current CHTC submission node, so Docker-based validation must be performed on another system with Docker installed.

## Option 1: Docker Validation (Recommended)

### Build and Test

```bash
# Build the Docker image
docker build -t redquanta .

# Test Module 1: Detect adaptive traits
docker run -v $(pwd)/results:/app/results redquanta \
    snakemake --configfile config/config_detect.yaml --cores 4

# Test Module 2: Evaluate performance
docker run -v $(pwd)/results:/app/results redquanta \
    snakemake evaluate_all --configfile config/config_evaluate.yaml --cores 2

# Check outputs
ls -la results/detect/
ls -la results/evaluate/
```

### Expected Outputs

**Module 1:**
- `results/detect/qst_results.csv` - Detection results for example traits

**Module 2:**
- `results/evaluate/*/tpr_fpr_matrix_*.csv` - TPR/FPR matrices
- `results/evaluate/combined_model_ranking.csv` - Model ranking

## Option 2: Local Conda Validation

If Docker is not available, validate using conda:

```bash
# Create and activate environment
conda env create -f environment.yml
conda activate redquanta

# Run dry-run first
snakemake --configfile config/config_detect.yaml --cores 4 -n

# Run Module 1
snakemake --configfile config/config_detect.yaml --cores 4

# Run Module 2 (reduced parameters)
snakemake evaluate_all --configfile config/config_evaluate.yaml --cores 2
```

## Option 3: HTCondor Validation

For validating HTCondor mode on CHTC or similar systems:

### Module 1: Single Trait Test

```bash
# Generate DAG for single trait
python htcondor/scripts/prepare_trait_dag.py \
    --trait-id L0MQ04 \
    --num-neutral 100 \
    --batch-size 50 \
    --sanity-check

# Submit DAG
condor_submit_dag results/dags/trait_L0MQ04.dag

# Monitor
condor_watch_q
```

### Module 2: Minimal Test

```bash
# Generate minimal evaluation DAG
python htcondor/scripts/prepare_perf_eval_dag.py \
    --chr autosomes \
    --num-repeats 100 \
    --num-neutral 100 \
    --batch-size 50 \
    --adaptive-qst 0.75,1.00 \
    --ve-ratios 1.0 \
    --output-dir results/perf_eval_test

# Submit
condor_submit_dag results/perf_eval_test/perf_eval_autosomes.dag
```

## Troubleshooting

### Common Issues

1. **Snakemake not found**: Ensure conda environment is activated
2. **R package errors**: Run `Rscript -e "library(abc)"` to verify
3. **Memory errors**: Reduce batch_size or num_sim in config
4. **Missing input files**: Check paths in config files

### Log Files

- Snakemake logs: `results/*/logs/*.log`
- HTCondor logs: `*.log`, `*.out`, `*.err` in output directories
