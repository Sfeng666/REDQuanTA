<div align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/logo/REDQuanTEA_logo_dark.png">
  <img alt="REDQuanTEA" src="docs/logo/REDQuanTEA_logo_light.png" width="100%">
</picture>

<h1>
<font color="#C41E3A">RED</font>Quan<font color="#C41E3A">TEA</font>
</h1>

<p><b>R</b>eplication-<b>E</b>nhanced <b>D</b>etection of <b>Quan</b>titative <b>T</b>raits <b>E</b>volving <b>A</b>daptively</p>
</div>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A workflow for detecting adaptive quantitative trait divergence (Q<sub>ST</sub>) from replicated measurements, using Approximate Bayesian Computation (ABC). It runs locally with Snakemake or on HTCondor (CHTC / OSPool).

Q<sub>ST</sub>–F<sub>ST</sub> comparisons lose power when extrinsic (non-genetic) variance deflates Q<sub>ST</sub>. REDQuanTEA uses biological replication plus ABC to separate V<sub>E</sub> from genetic variance, then compares each trait to a trait-specific dynamic neutral Q<sub>ST</sub> threshold matched on sample structure, estimator, chromosome, and V<sub>E</sub>/V<sub>G</sub>.

## Workflow overview

![QST Workflow](data/reference/Figure%201_QST_workflow_both_modes.png)

### Detection Module

1. Calculate variance components (method of moments) with the ridge floor.
2. Estimate trait Q<sub>ST</sub> with loclinear ABC.
3. Build a neutral Q<sub>ST</sub> distribution from empirical F<sub>ST</sub> values at the same V<sub>E</sub>/V<sub>G</sub>.
4. Compare trait Q<sub>ST</sub> to the 95th percentile of that null.
5. Call a trait adaptive if Q<sub>ST</sub> exceeds the threshold.

### Design Module

1. Simulate traits with known adaptive Q<sub>ST</sub> and V<sub>E</sub>/V<sub>G</sub>.
2. Compute true positive rate (TPR) and false positive rate (FPR).
3. Rank summary-statistic combinations.
4. Compare power across sample structures (optional).

Manuscript figure sources are listed in [docs/MANUSCRIPT_FIGURES.md](docs/MANUSCRIPT_FIGURES.md).

## Installation

```bash
git clone https://github.com/sfeng666/REDQuanTEA.git
cd REDQuanTEA

conda env create -f environment.yml
conda activate redquantea
```

On CHTC submit nodes, load Snakemake from the [htc-modules](https://chtc.cs.wisc.edu/uw-research-computing/htc-modules) collection if conda is not already set up.

Pack an R environment for HTCondor workers once:

```bash
bash htcondor/scripts/pack_r_env.sh
bash htcondor/scripts/sync_r_env_to_staging.sh
```

## Quick start (local / Snakemake)

Short check (2 traits, fewer simulations):

```bash
bash scripts/validate_local.sh
# or:
snakemake --configfile config/config_detect_smoke.yaml --cores 4
snakemake evaluate_all --configfile config/config_evaluate_smoke.yaml --cores 4
```

Detection Module on the full example table:

```bash
snakemake --configfile config/config_detect.yaml --cores 4 -n
snakemake --configfile config/config_detect.yaml --cores 4
# Output: results/detect/qst_results.csv
```

Design Module (reduced parameters for a workstation):

```bash
snakemake evaluate_all --configfile config/config_evaluate.yaml --cores 4
# Output: results/evaluate/combined_model_ranking.csv
```

## Quick start (HTCondor / CHTC)

Detection Module, fused one-job-per-trait DAG:

```bash
TRAIT_VALUES=data/example/trait_values_smoke.csv \
RESULTS_DIR=results/detect_chtc \
OUTPUT_DAG=results/dags/detect.dag \
NUM_NEUTRAL=100 NUM_SIM=10000 MAX_TRAITS=2 DRY_RUN=1 \
bash htcondor/scripts/submit_detection.sh
```

Set `DRY_RUN=0` (default) to submit. After jobs finish:

```bash
bash htcondor/scripts/aggregate_detection.sh results/detect_chtc
```

Design Module, one summary-stat combination:

```bash
printf 'QST\tratioVbetweenVtotal\n' > /tmp/combo.txt
COMBO_FILE=/tmp/combo.txt \
OUTPUT_DIR=results/design_chtc \
NUM_REPEATS=100 NUM_NEUTRAL=100 BATCH_SIZE=50 DRY_RUN=1 \
bash htcondor/scripts/submit_design.sh
```

After jobs finish, aggregate TPR/FPR matrices:

```bash
bash workflow/scripts/run_aggregate_perf_eval_multicombo_fast.sh results/design_chtc
```

## Input files

| File | Description |
|------|-------------|
| `sample_structure.csv` | Population / strain / replicate layout |
| `trait_values.csv` | Trait measurements per strain / replicate |
| `qst_neutral_autosomes.txt` | Neutral F<sub>ST</sub> values (autosomes) |
| `qst_neutral_chrX.txt` | Neutral F<sub>ST</sub> values (X chromosome) |

Example copies live in `data/example/`. Sample structure:

```csv
population,strain,replicate
pop1,strain1,rep1
pop1,strain1,rep2
pop1,strain2,rep1
```

Trait values (wide format: one row per trait, sample columns after `chr`):

```csv
trait_id,chr,1,2,3
trait001,autosomes,0.523,0.541,0.510
```

## Output files

| File | Description |
|------|-------------|
| `qst_results.csv` | Detection results (`trait_id`, `chr`, `QST`, `adaptive`) |
| `tpr_fpr_matrix_*.csv` | TPR/FPR matrices per chromosome |
| `combined_model_ranking.csv` | Model performance ranking |

Blank `QST` / `adaptive` fields are not always a pipeline error. If observed ANOVA components are uninformative (for example both genetic variances non-positive), ABC can return NA. Threshold columns can still be written; treat the blank estimate as "uninformative under this sample/trait pattern".

## Extra: compare several trait types

Not required by the manuscript, but useful if you run Detection on more than one collection. The scripts take any set of directories that contain `summary/qst_results_abc.csv`:

```bash
Rscript workflow/scripts/plot_comparison_across_trait_types.R \
  --results-root results/detect_collections \
  --output-dir results/detect_collections/comparison

Rscript workflow/scripts/compare_adaptive_fst_vs_qst.R \
  --results-root results/detect_collections \
  --output-dir results/detect_collections/comparison \
  --fst-autosomes data/example/qst_neutral_autosomes.txt \
  --fst-chrx data/example/qst_neutral_chrX.txt
```

Optional `--trait-types a,b,c` and `--labels-tsv` (columns `trait_type`, `label`) override auto-discovery and plot labels.

## Sample-structure comparison (optional)

When `sample_structures` is set in `config_evaluate.yaml`, Snakemake evaluates each structure. After HTCondor Design runs that write `n2_i{ind}_r{rep}/` trees:

```bash
bash workflow/scripts/run_aggregate_sample_structure_perf_eval.sh \
  path/to/sample_struct_results \
  0.8 0.1,1.0,10.0 0.95

bash workflow/scripts/run_plot_sample_structure_comparison.sh \
  path/to/sample_struct_results \
  "QST,ratioVbetweenVtotal"
```

## Documentation

- [README_details.md](README_details.md): parameters, HTCondor setup, troubleshooting
- [docs/MANUSCRIPT_FIGURES.md](docs/MANUSCRIPT_FIGURES.md): manuscript figures to generating scripts
- [VALIDATION.md](VALIDATION.md): local and CHTC checks
- [HTCondor documentation](https://htcondor.readthedocs.io/)

## Citation

If you use REDQuanTEA, please cite:

> Feng, S., & Pool, J. E. (in preparation). *REDQuanTEA: Replication-Enhanced Detection of Quantitative Traits Evolving Adaptively*. Laboratory of Genetics, University of Wisconsin-Madison.

Until the paper is published, you may also cite this repository:

> Feng, S., & Pool, J. E. (2026). *REDQuanTEA* [Software]. GitHub. https://github.com/sfeng666/REDQuanTEA

## License

MIT. See [LICENSE](LICENSE).
