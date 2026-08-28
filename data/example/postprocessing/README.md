# Post-processing examples (Git-friendly)

Small fixtures for optional local steps documented in [README.md](../../README.md) and [README_details.md](../../README_details.md). **Total tree size ~200 KB** (12 files of RData + CSV/PDF examples).

## Layout

| Step | Script | Example input | Example output |
|------|--------|---------------|----------------|
| Sample-structure aggregation | [workflow/scripts/run_aggregate_sample_structure_perf_eval.sh](../../workflow/scripts/run_aggregate_sample_structure_perf_eval.sh) | [sample_struct_aggregation/input/](sample_struct_aggregation/input/) | [sample_struct_aggregation/output/](sample_struct_aggregation/output/) |
| Sample-structure plots | [workflow/scripts/run_plot_sample_structure_comparison.sh](../../workflow/scripts/run_plot_sample_structure_comparison.sh) | [sample_struct_plots/input/](sample_struct_plots/input/) | [sample_struct_plots/output/](sample_struct_plots/output/) |
| Multi-combo perf-eval aggregation | [workflow/scripts/run_aggregate_perf_eval_multicombo_fast.sh](../../workflow/scripts/run_aggregate_perf_eval_multicombo_fast.sh) | HTCondor `autosomes/` + `chrX/` batch `.RData` (see below) | [perf_eval_multicombo/output/](perf_eval_multicombo/output/) |
| Publication model ranking | [workflow/scripts/run_perf_eval_publication_ranking.sh](../../workflow/scripts/run_perf_eval_publication_ranking.sh) | [perf_eval_ranking/input/combined_model_ranking.csv](perf_eval_ranking/input/combined_model_ranking.csv) | [perf_eval_ranking/output/](perf_eval_ranking/output/) |

## File inventory (bytes, approximate)

### `sample_struct_aggregation/` (~130 KB)

**Input** — one structure `n2_i4_r3`, subset of a real validation run (1 batch file per condition):

- `input/n2_i4_r3/autosomes/neutral_ratio_{1,2,3}/neutral_batch_1.RData` — 3 × ~10 KB
- `input/n2_i4_r3/autosomes/adaptive_q0_80_r{1,2,3}/adaptive_batch_1.RData` — 3 × ~8 KB
- `input/n2_i4_r3/chrX/neutral_ratio_1/neutral_batch_1.RData` — ~10 KB
- `input/n2_i4_r3/chrX/adaptive_q0_80_r1/adaptive_batch_1.RData` — ~8 KB

**Output** — aggregated matrices (regenerate with the aggregation script above):

- `output/tpr_fpr_matrix_autosomes.csv` — 161 B
- `output/tpr_fpr_matrix_chrX.csv` — 120 B

### `sample_struct_plots/` (~56 KB)

**Input** — two structures with completed `tpr_fpr_matrix_*.csv` (plotting does not read `.RData`):

- `input/n2_i4_r3/{autosomes,chrX}/tpr_fpr_matrix_*.csv` — 4 files, ~160 B each
- `input/n2_i6_r3/{autosomes,chrX}/tpr_fpr_matrix_*.csv` — 4 files, ~160 B each

**Output** — regenerate with the plotting script:

- `output/power_comparison_combined_QST.ratioVbetweenVtotal.pdf` — ~6.4 KB
- `output/power_comparison_combined_per_pop_QST.ratioVbetweenVtotal.pdf` — ~6.4 KB
- `output/power_summary_table_{autosomes,chrX}_QST.ratioVbetweenVtotal.txt` — ~125 B each

### `perf_eval_multicombo/` (~9 KB)

Illustrates **output file formats** only (full runs produce thousands of `.RData` files and are not checked in).

- `combinations_smoke_3.txt` — 75 B (3 tab-separated combos)
- `output/tpr_fpr_matrix_autosomes_combo_0001.csv` — 314 B (schema reference)
- `output/combined_model_ranking.csv` — 253 B (schema reference)

### `perf_eval_ranking/` (~2 KB)

Publication-friendly ranking from `combined_model_ranking.csv` (internal stat names → manuscript labels; drops skew/kurtosis combos).

- `input/combined_model_ranking.csv` — top 5 rows from a 2047-combo validation run
- `output/combined_model_ranking_publication.csv` — same rows with publication stat names
- `output/Table_model_ranking.txt` — two-column table at V_E/V_G = 1.0
- `output/Table_model_ranking_legend.txt` — table caption text

Regenerate the example outputs:

```bash
Rscript workflow/scripts/format_model_ranking.R \
  data/example/postprocessing/perf_eval_ranking/input 1.0
# then move combined_model_ranking_publication.csv and Table_* into output/
```

Real multi-combo input layout under a perf-eval root:

```text
perf_eval_root/
  combinations_all_nonempty_11stats.txt
  autosomes/neutral_ratio_1/neutral_batch_*.RData
  autosomes/adaptive_q0_50_r1/adaptive_batch_*.RData
  chrX/...
```

Each batch `.RData` has `result$is_multi_combo == TRUE` and a `combo` column when using the fast all-combo worker.

## Reproduce examples

From the REDQuanTEA repository root with `redquantea` (or any env with `ggplot2`, `dplyr`, `tidyr`, `viridis`):

```bash
# Aggregation example
bash workflow/scripts/run_aggregate_sample_structure_perf_eval.sh \
  data/example/postprocessing/sample_struct_aggregation/input

# Plots example
bash workflow/scripts/run_plot_sample_structure_comparison.sh \
  data/example/postprocessing/sample_struct_plots/input \
  "QST,ratioVbetweenVtotal" \
  data/example/postprocessing/sample_struct_plots/output
```

## Full validation results (not in git)

Production-scale sample-structure comparison (12 structures, both chromosomes):

`htcondor/results/validation_sample_struct_QST_ratioVbetweenVtotal/`

After HTCondor completes, run aggregation then plots with the same scripts and `adaptive_qst=0.8`, `ve_ratios=0.1,1.0,10.0`.

Fast 2047-combo perf-eval (autosomes + chrX):

`htcondor/results/validation_perf_fast_2047/`

```bash
# If per-combo CSVs already exist (e.g. after marula aggregation):
bash workflow/scripts/run_perf_eval_publication_ranking.sh \
  htcondor/results/validation_perf_fast_2047
```
