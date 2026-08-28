# REDQuanTEA validation

## Local (Snakemake)

Short smoke (2 traits, 2000 ABC simulations, tiny Design grid):

```bash
conda activate redquantea
bash scripts/validate_local.sh
```

Expected Detection output: `results/detect_smoke/qst_results.csv`.
Expected Design output: `results/evaluate_smoke/combined_model_ranking.csv` and `tpr_fpr_matrix_*.csv`.

Fuller local configs (still reduced relative to CHTC):

```bash
snakemake --configfile config/config_detect.yaml --cores 4
snakemake evaluate_all --configfile config/config_evaluate.yaml --cores 4
```

`config/config_detect.yaml` uses `data/example/trait_values.csv`, which is the full example protein table. Use the smoke config on a laptop.

## HTCondor smoke

Generate only (`DRY_RUN=1`):

```bash
TRAIT_VALUES=data/example/trait_values_smoke.csv \
RESULTS_DIR=results/detect_chtc_smoke \
OUTPUT_DAG=results/dags/detect_smoke.dag \
NUM_NEUTRAL=20 NUM_SIM=2000 MAX_TRAITS=2 BATCH_SIZE=20 DRY_RUN=1 \
bash htcondor/scripts/submit_detection.sh

printf 'QST\tratioVbetweenVtotal\n' > /tmp/combo.txt
COMBO_FILE=/tmp/combo.txt \
OUTPUT_DIR=results/design_chtc_smoke \
NUM_REPEATS=20 NUM_NEUTRAL=20 BATCH_SIZE=20 NUM_SIM=2000 \
CHR=autosomes DRY_RUN=1 \
bash htcondor/scripts/submit_design.sh
```

Set `DRY_RUN=0` to submit. After Detection jobs finish:

```bash
bash htcondor/scripts/aggregate_detection.sh results/detect_chtc_smoke
```

After Design jobs finish, the DAG POST script writes `tpr_fpr_matrix_*.csv`. You can re-run:

```bash
bash workflow/scripts/run_aggregate_perf_eval_multicombo_fast.sh results/design_chtc_smoke
```

## Docker

```bash
docker build -t redquantea .
docker run -v $(pwd)/results:/app/results redquantea \
  snakemake --configfile config/config_detect_smoke.yaml --cores 2
```
