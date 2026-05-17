# REDQuanTA Detailed Documentation

This document provides comprehensive documentation for REDQuanTA, including full parameter descriptions, HTCondor setup instructions, and troubleshooting guides.

## Table of Contents

1. [Installation](#installation)
2. [Configuration Files](#configuration-files)
3. [Module 1: Detect Adaptive Traits](#module-1-detect-adaptive-traits)
4. [Module 2: Evaluate Performance](#module-2-evaluate-performance)
5. [HTCondor Setup](#htcondor-setup)
6. [FST Input Options](#fst-input-options)
7. [Performance Benchmarks](#performance-benchmarks)
8. [Troubleshooting](#troubleshooting)
9. [Technical Details](#technical-details)

---

## Installation

### Prerequisites

- Conda/Mamba (recommended: [Miniforge](https://github.com/conda-forge/miniforge))
- For HTCondor execution: Access to an HTCondor pool

### Environment Setup

```bash
# Using conda
conda env create -f environment.yml
conda activate redquanta

# Using mamba (faster)
mamba env create -f environment.yml
mamba activate redquanta
```

### Verify Installation

```bash
# Check Snakemake
snakemake --version

# Check R packages
Rscript -e "library(abc); library(ggplot2); cat('R packages OK\n')"
```

---

## Configuration Files

### config/config_detect.yaml (Module 1)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sample_structure` | `data/example/sample_structure.csv` | Path to sample structure file |
| `trait_values` | `data/example/trait_values.csv` | Path to trait values file |
| `fst_input_mode` | `direct` | FST input mode: `direct`, `from_vcf`, or `from_simulation` |
| `fst_autosomes` | `data/example/qst_neutral_autosomes.txt` | Path to autosome FST values |
| `fst_chrX` | `data/example/qst_neutral_chrX.txt` | Path to X chromosome FST values |
| `num_neutral` | 100 | Number of neutral FST values (reduce for testing) |
| `num_sim` | 100000 | ABC simulations per estimation |
| `batch_size` | 50 | FST values per batch job |
| `threshold_percentile` | 0.95 | Threshold for adaptive detection |
| `summary_stats` | `QST,ratioVbetweenVtotal` | Summary statistics for ABC |
| `chromosomes` | `[autosomes, chrX]` | Chromosome types to analyze |
| `output_dir` | `results/detect` | Output directory |

### config/config_evaluate.yaml (Module 2 - Local)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `adaptive_qst_levels` | `[0.50, 0.75, 1.00]` | Adaptive Q<sub>ST</sub> levels to test |
| `ve_ratios` | `[0.1, 1.0, 10.0]` | V<sub>E</sub>/V<sub>G</sub> ratios |
| `num_repeats` | 100 | Repeats per condition (reduced for local) |
| `num_neutral` | 100 | Neutral simulations (reduced for local) |
| `num_sim` | 100000 | ABC simulations per estimation |
| `batch_size` | 50 | Values per batch |
| `summary_stats_combos` | `[QST,ratioVbetweenVtotal, ...]` | Summary stat combinations to evaluate |

### config/config_evaluate_full.yaml (Module 2 - HTCondor)

Full-scale parameters for HTCondor execution:

| Parameter | Full Value | Description |
|-----------|------------|-------------|
| `adaptive_qst_levels` | 0.50–1.00 (0.05 steps) | 11 levels |
| `ve_ratios` | `[0.01, 0.1, 1.0, 10.0, 100.0]` | 5 ratios |
| `num_repeats` | 10000 | Full evaluation |
| `num_neutral` | 10000 | Full neutral distribution |
| `batch_size` | 1000 | Optimal for HTCondor |

---

## Module 1: Detect Adaptive Traits

### Workflow Steps

1. **prepare_obs_stats**: Calculate observed variance components per trait
2. **estimate_trait_qst**: ABC estimation of trait Q<sub>ST</sub>
3. **estimate_neutral_qst**: ABC estimation for neutral FST values (batched)
4. **aggregate_trait**: Compare trait Q<sub>ST</sub> to neutral threshold
5. **aggregate_all_traits**: Combine results into final CSV

### Running with Snakemake

```bash
# Dry run
snakemake --configfile config/config_detect.yaml --cores 4 -n

# Full run
snakemake --configfile config/config_detect.yaml --cores 4

# With custom parameters
snakemake --configfile config/config_detect.yaml \
    --config num_neutral=1000 batch_size=100 \
    --cores 8
```

### Output: qst_results.csv

| Column | Description |
|--------|-------------|
| trait_id | Trait identifier |
| chr | Chromosome type |
| QST | Estimated trait Q<sub>ST</sub> |
| threshold_percentile | Percentile used (0.95) |
| threshold_value | Neutral Q<sub>ST</sub> at threshold |
| adaptive | "yes" if Q<sub>ST</sub> > threshold |

---

## Module 2: Evaluate Performance

### Workflow Steps

1. **evaluate_neutral_qst**: Estimate Q<sub>ST</sub> for neutral FST values
2. **evaluate_adaptive_qst**: Estimate Q<sub>ST</sub> for simulated adaptive traits
3. **aggregate_perf_eval**: Calculate TPR/FPR matrices
4. **combined_model_ranking**: Rank summary stat combinations
5. **plot_tpr_performance**: Generate TPR plots

### Running with Snakemake

```bash
# Local (reduced parameters)
snakemake evaluate_all --configfile config/config_evaluate.yaml --cores 4

# Full evaluation (HTCondor recommended)
snakemake evaluate_all --configfile config/config_evaluate_full.yaml --cores 16
```

### Summary Statistics Combinations

| Combination | Description |
|-------------|-------------|
| `QST,ratioVbetweenVtotal` | Q<sub>ST</sub> + V<sub>among</sub>/(V<sub>among</sub>+V<sub>within</sub>+V<sub>E</sub>) (default benchmark combo) |
| `QST,ratioVext` | Q<sub>ST</sub> + V<sub>E</sub>/V<sub>G</sub> ratio |
| `QST,ext_sd` | Q<sub>ST</sub> + extrinsic SD |
| `QST` | Q<sub>ST</sub> only |

---

## HTCondor Setup

### Prerequisites

1. Access to an HTCondor pool (e.g., UW-Madison CHTC, OSG)
2. Submit node with HTCondor client installed

### R Environment Setup

HTCondor worker jobs expect **`htcondor/env/r_env.tar.gz`** (a **conda-pack** of the conda env defined in **`htcondor/env/r_qst.yml`**: `r-base`, `r-e1071`, `r-abc`). Local Snakemake can still use the full root **`environment.yml`** for plotting and orchestration.

From the repository root, build the tarball on a machine with enough disk quota:

```bash
bash htcondor/scripts/pack_r_env.sh
```

The output is `htcondor/env/r_env.tar.gz` (listed in `.gitignore` because of size; add it to your release artifacts or transfer out-of-band as needed).

### Submitting Jobs

#### Module 1: Single Trait

```bash
python htcondor/scripts/prepare_trait_dag.py \
    --trait-id L0MQ04 \
    --num-neutral 10000 \
    --batch-size 1000

condor_submit_dag results/dags/trait_L0MQ04.dag
```

#### Module 1: Multiple Traits

```bash
python htcondor/scripts/generate_all_dags.py \
    --max-traits 100 \
    --num-neutral 10000 \
    --batch-size 1000

condor_submit_dag results/dags/all_traits.dag
```

#### Module 2: Performance Evaluation

```bash
python htcondor/scripts/prepare_perf_eval_dag.py \
    --chr both \
    --num-repeats 10000 \
    --output-dir results/perf_eval

condor_submit_dag results/perf_eval/perf_eval_autosomes.dag
condor_submit_dag results/perf_eval/perf_eval_chrX.dag
```

### Optional Post-Processing (Local)

These steps are **not** run automatically by HTCondor DAGs. Run them on the
submit node or your workstation after the corresponding jobs finish.

#### Sample-structure aggregation + power plots

HTCondor Module 2 DAGs write batch `.RData` under `neutral_ratio_{i}/` and
`adaptive_q{qst}_r{i}/`. The DAG POST script calls `aggregate_perf_eval.R`, which
must read those **directories** (not legacy flat log basenames). If matrices are
all `NA`, re-run aggregation locally with the script below.

**Step A — aggregate** (same `adaptive_qst` / `ve_ratios` as in your DAG):

```bash
bash workflow/scripts/run_aggregate_sample_structure_perf_eval.sh \
  /path/to/sample_struct_results \
  0.8 0.1,1.0,10.0 0.95
```

**Step B — plot** (requires numeric `tpr_fpr_matrix_*.csv` from step A):

```bash
bash workflow/scripts/run_plot_sample_structure_comparison.sh \
  /path/to/sample_struct_results \
  "QST,ratioVbetweenVtotal" \
  /path/to/sample_struct_results/plots
```

Example (CHTC validation tree):

```bash
bash workflow/scripts/run_aggregate_sample_structure_perf_eval.sh \
  htcondor/results/validation_sample_struct_QST_ratioVbetweenVtotal
bash workflow/scripts/run_plot_sample_structure_comparison.sh \
  htcondor/results/validation_sample_struct_QST_ratioVbetweenVtotal
```

**Outputs**

- Per structure/chromosome: `tpr_fpr_matrix_{autosomes,chrX}.csv`, `*_heatmap.pdf`
- Plots dir: `power_comparison_combined_*.pdf`, `power_summary_table_*.txt`

**Examples (checked in, ~200 KB total):** [data/example/postprocessing/README.md](../data/example/postprocessing/README.md)

| Path | Role |
|------|------|
| [sample_struct_aggregation/input/](data/example/postprocessing/sample_struct_aggregation/input/) | Mini `n2_i4_r3` batch `.RData` tree |
| [sample_struct_aggregation/output/](data/example/postprocessing/sample_struct_aggregation/output/) | Golden `tpr_fpr_matrix_*.csv` |
| [sample_struct_plots/input/](data/example/postprocessing/sample_struct_plots/input/) | Two structures × 2 chr matrices |
| [sample_struct_plots/output/](data/example/postprocessing/sample_struct_plots/output/) | Golden PDFs + TPR tables |

The R script `workflow/scripts/plot_sample_structure_comparison.R` auto-discovers
`n2_i*_r*` subdirectories.

#### Fast multi-combo perf-eval aggregation

Use when each HTCondor job wrote **many** summary-stat combinations into batch
`.RData` files (e.g. `prepare_perf_eval_all_combos_fast_dag.py` with a
`combinations*.txt` file). This is much faster than calling
`aggregate_perf_eval.R` once per combo.

```bash
bash workflow/scripts/run_aggregate_perf_eval_multicombo_fast.sh \
  /path/to/perf_eval_root \
  /path/to/perf_eval_root/combinations_all_nonempty_11stats.txt
```

**Inputs:** `autosomes/` and `chrX/` with `neutral_ratio_*` and `adaptive_q*_*`
subdirectories containing `neutral_batch_*.RData` / `adaptive_batch_*.RData`.

**Outputs:**

- Per combo: `tpr_fpr_matrix_{chr}_combo_XXXX.csv` under each chromosome dir
- `combined_model_ranking.csv` at the results root (via `generate_combined_ranking.R`)

#### Publication model ranking

After per-combo matrices exist (and optionally `combined_model_ranking.csv`):

```bash
bash workflow/scripts/run_perf_eval_publication_ranking.sh /path/to/perf_eval_root
# optional second argument: V_E/V_G for Table_model_ranking.txt (default 1.0)
```

Or set `RUN_PUBLICATION_RANKING=1` when calling
`run_aggregate_perf_eval_multicombo_fast.sh` to chain aggregation and publication
formatting. Use `SKIP_COMBINED=1` if `combined_model_ranking.csv` already exists.

**Outputs:**

- `combined_model_ranking_publication.csv` — full ranking with publication stat labels
- `Table_model_ranking.txt` — tab-delimited two-column table (TPR at one V_E/V_G)
- `Table_model_ranking_legend.txt` — caption / abbreviation key

R scripts: `generate_combined_ranking.R` (reads `model` column from per-combo CSVs),
`format_model_ranking.R` (rename stats, filter skew/kurtosis).

**Example (5 rows, ~2 KB):** [data/example/postprocessing/perf_eval_ranking/](../data/example/postprocessing/perf_eval_ranking/)

**Environment overrides:** `ADAPTIVE_QST`, `VE_RATIOS`, `THRESHOLD`, `RSCRIPT`
(see script header). Defaults match the standard 11-level adaptive Q<sub>ST</sub>
grid and five V<sub>E</sub>/V<sub>G</sub> ratios used in full perf-eval DAGs.

For single-combo Module 2 runs, the DAG POST script still uses
`aggregate_perf_eval.R`; use the fast script only when many combos were scored
per batch job.

**Example output formats (synthetic, ~9 KB):** [data/example/postprocessing/perf_eval_multicombo/](../data/example/postprocessing/perf_eval_multicombo/)

### HTCondor Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--batch-size` | 1000 | FST values per job (optimal for 5-min policy) |
| `--num-sim` | 100000 | ABC simulations |
| `--priority` | None | Job priority (higher = first) |
| `--sanity-check` | False | Enable diagnostic output |

### Monitoring Jobs

```bash
condor_q                    # List running jobs
condor_q -hold              # List held jobs
condor_q -dag               # DAG status
condor_watch_q              # Real-time monitoring
tail -f results/dags/*.dagman.out  # DAG progress
```

---

## FST Input Options

### Option 1: Direct (default)

Provide pre-computed FST files:

```yaml
fst_input_mode: "direct"
fst_autosomes: "data/example/qst_neutral_autosomes.txt"
fst_chrX: "data/example/qst_neutral_chrX.txt"
```

### Option 2: From VCF

Generate FST from VCF file:

```yaml
fst_input_mode: "from_vcf"
vcf_file: "path/to/variants.vcf.gz"
vcf_populations: "pop1,pop2"
```

### Option 3: From Simulation

Generate FST using ms coalescent simulation:

```yaml
fst_input_mode: "from_simulation"
ms_config: "path/to/ms_config.yaml"
```

---

## Performance Benchmarks

### Batch Size Optimization

Based on HTCondor benchmarking:

| Batch Size | Avg Job Time | Time per FST | Recommendation |
|------------|--------------|--------------|----------------|
| 100 | 13.6 min | 8.20s | Too short |
| 500 | 69.0 min | 8.28s | Acceptable |
| **1000** | **109.7 min** | **6.58s** | **Optimal** |
| 2000 | 238.1 min | 7.14s | Too long |

**Recommendation**: Use `batch_size: 1000` for HTCondor, `batch_size: 50` for local testing.

### Memory Usage

- Per-FST memory: ~100 MB peak
- Optimized to release memory between FST estimations
- HTCondor jobs request 4GB, auto-retry on memory exceeded

---

## Troubleshooting

### Common Issues

#### 1. Snakemake: Missing Input Files

```
Error: Missing input files for rule prepare_obs_stats
```

**Solution**: Verify input files exist and paths in config are correct.

#### 2. HTCondor: Jobs Held

```bash
condor_q -hold
# Check hold reason
condor_q -af HoldReason
```

Common reasons:
- Memory exceeded: Jobs auto-retry on different nodes
- Transfer failures: Auto-retry up to 5 times

#### 3. R Package Loading Errors

```
Error in library(abc): there is no package called 'abc'
```

**Solution**: Reinstall environment:
```bash
conda env remove -n redquanta
conda env create -f environment.yml
```

#### 4. Empty Results

Check logs:
```bash
# Snakemake
cat results/detect/logs/*.log

# HTCondor
cat results/dags/*.err
```

### Getting Help

1. Check log files for specific error messages
2. Run with `--sanity-check` flag for diagnostic output
3. For Snakemake: Use `-p` flag to print commands

---

## Technical Details

### ABC Estimation Parameters

- **Simulations**: 100,000 per estimation
- **Method**: Local linear (`loclinear`, default) or `rejection` use the vectorized estimator in `abc_fast_estimators.R` (numerically equivalent to `abc::abc` on a fixed pool). `neuralnet` and `ridge` still use **`abc::abc`**. Override with `QST_ABC_METHOD`.
- **Tolerance**: Dynamically adjusted for sufficient accepted samples
- **Summary statistics**: Configurable; default: QST, ratioVbetweenVtotal

### Variance Components

Using Method of Moments (MoM):

- V<sub>among</sub>: Among-population variance
- V<sub>within</sub>: Within-population, among-strain variance
- V<sub>E</sub>: Residual (environmental) variance

Q<sub>ST</sub> = V<sub>among</sub> / (V<sub>among</sub> + 2 × V<sub>within</sub>)

### Memory Optimization

The workflow uses vectorized operations and explicit garbage collection to minimize memory usage:

- Pre-computed indices for variance calculations
- `gc(reset=TRUE)` after each FST to return memory to OS
- Reduces per-FST memory growth from ~36 MB to ~0.3 MB

---

## References

- HTCondor: https://htcondor.readthedocs.io/
- Snakemake: https://snakemake.readthedocs.io/
- ABC package: https://cran.r-project.org/package=abc
