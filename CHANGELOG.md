# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Ridge variance floor (`FLOOR_POLICY=ridge_floor`, `FLOOR_ALPHA=0.1`) as the Detection and Design Module default, matching the manuscript.
- Detection Module fused HTCondor DAG (`prepare_detection_dag.py`, `submit_detection.sh`) and Design Module DAG (`prepare_design_dag.py`, `submit_design.sh`).
- Local smoke configs: `config/config_detect_smoke.yaml`, `config/config_evaluate_smoke.yaml`, `data/example/trait_values_smoke.csv`.
- Trait-type comparison scripts that take `--results-root` rather than a fixed collection list: `compare_adaptive_fst_vs_qst.R`, `plot_comparison_across_trait_types.R`.
- Manuscript figure map in `docs/MANUSCRIPT_FIGURES.md` and dark/light logos under `docs/logo/`.
- Publication perf-eval ranking: `format_model_ranking.R` (synced from qst_workflow), `run_perf_eval_publication_ranking.sh`, and example fixtures under `data/example/postprocessing/perf_eval_ranking/`.
- `generate_combined_ranking.R`: use internal `model` column from per-combo CSVs (not `combo_*` filenames) when building `combined_model_ranking.csv`.
- Optional post-processing wrappers: `run_aggregate_sample_structure_perf_eval.sh`, `run_plot_sample_structure_comparison.sh`, `run_aggregate_perf_eval_multicombo_fast.sh`; examples under `data/example/postprocessing/` (~200 KB).
- `workflow/scripts/aggregate_perf_eval_multicombo_fast.R` for local multi-combo perf-eval aggregation.
- Optional script `workflow/scripts/perf_eval_all_combos_fast_cell_v2.R` for large all-subset performance evaluation (same sourced ABC core as Design Module; not part of the default Snakemake graph).
- `workflow/scripts/abc_fast_estimators.R`: shared vectorized loclinear / rejection ABC (same outputs as `abc::abc` on a fixed simulation pool). Loaded from `qst_abc_sim.R`; used automatically when `QST_ABC_METHOD` is `loclinear` or `rejection`.
- Optional script `workflow/scripts/benchmark_loclinear_fast_vs_abc.R`: on a shared simulation pool, compares `fast_abc_estimate_multi` to `abc::abc` with `method = "loclinear"` (numeric agreement and relative timings).
- `prepare_obs_stats.R` optional sixth argument: path to write observed V_E/V_G (`ratioVext`) as text for neutral-batch ABC jobs.
- Output `{trait}_ratioVext.txt` from Detection Module prep when using Snakemake.
- `htcondor/env/r_qst.yml` and `htcondor/env/r_env.tar.gz`: conda recipe and conda-packed R bundle for HTCondor workers; `htcondor/scripts/pack_r_env.sh` rebuilds `r_env.tar.gz` from `r_qst.yml`. Tarball remains gitignored by default due to size.

### Changed

- Project name is REDQuanTEA (Replication-Enhanced Detection of Quantitative Traits Evolving Adaptively). User-facing text uses Detection Module and Design Module.
- Canonical ABC engine lives in `workflow/scripts/qst_abc_sim.R` (no `_opt` suffix). Default method remains loclinear.
- Design Module Snakemake jobs write `neutral_ratio_{i}/` and `adaptive_q{qst}_r{i}/` directories so `aggregate_perf_eval.R` can find them.
- **ABC core (`workflow/scripts/qst_abc_sim.R`)** aligned with the soft-clip / prior-flooring reference implementation: expanded summary-stat pool (including ratio summaries), soft-clipped negative ANOVA variance components in simulation, **prior_floor** on genetic SD priors (0.1× sum of observed SDs), `QST_ABC_METHOD` and `QST_ABC_TOL_NUMERATOR` environment controls, and `run_abc_qst` / `qst_mean_from_abc` wiring consistent with `abc::abc`.
- **Vectorized ABC**: when `QST_ABC_METHOD` is `loclinear` or `rejection`, `run_abc_qst` (all single-combo paths) and `estimate_qst_abc_multi` use `abc_fast_estimators.R` instead of `abc::abc`; `neuralnet` and `ridge` still use `abc::abc`.
- **Default summary-stat combination** for Detection Module and Design Module is now **`QST,ratioVbetweenVtotal`** (configs, Snakemake rule defaults, HTCondor wrappers, and CLI fallback in `qst_abc_sim.R`).
- **Default `QST_ABC_METHOD`** in `qst_abc_sim.R` is now **`loclinear`** (was `neuralnet`); set the environment variable to override.
- **`prepare_obs_stats.R`**: observed summaries now include `ratioVbetweenVext`, `ratioVbetweenVtotal`, `ratioVwithinVtotal`, and `ratioVextVtotal`, matching the benchmark pipeline used for performance evaluation.
- **`README_details.md`**: default combo documentation updated to match.

### Fixed

- **`aggregate_perf_eval.R`**: load batch `.RData` from `neutral_ratio_{i}/` and `adaptive_q{qst}_r{i}/` directories (fixes all-`NA` `tpr_fpr_matrix_*.csv` after HTCondor Design Module).
- **Detection Module neutral batches (`batch_neutral`)**: Snakemake now passes **ratioVext** (V_E/V_G from the trait), not environmental SD (`ext_sd`), as the third argument—matching `generate_neutral_obs_stats` semantics.
- **HTCondor `prepare_trait_dag.py`**: local prep calls `prepare_obs_stats.R` with **`sample_structure` before `trait_values`**; DAG generation uses `SCRIPT_DIR` (not an undefined `SCRIPTS_DIR`) for `CONFIG` and `.sub` paths.
- **HTCondor `prepare_perf_eval_dag.py`**: DAG generation uses `SCRIPT_DIR` for `CONFIG` and `.sub` paths (same `SCRIPTS_DIR` NameError class of bug).
