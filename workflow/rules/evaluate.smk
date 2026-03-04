# REDQuanTA Module 2: Performance Evaluation Rules
#
# Workflow:
#   1. For each (adaptive_qst, ve_ratio, chr, summary_stats) combination:
#      - Estimate neutral QST distribution
#      - Estimate adaptive QST
#      - Calculate TPR/FPR
#   2. Aggregate into TPR/FPR matrices
#   3. Generate combined model ranking

import os
from pathlib import Path

# Get config values
ADAPTIVE_QST_LEVELS = config.get("adaptive_qst_levels", [0.50, 0.75, 1.00])
VE_RATIOS = config.get("ve_ratios", [0.1, 1.0, 10.0])
NUM_REPEATS = config.get("num_repeats", 100)
NUM_NEUTRAL = config.get("num_neutral", 100)
NUM_SIM = config.get("num_sim", 100000)
BATCH_SIZE = config.get("batch_size", 50)
THRESHOLD_PERCENTILE = config.get("threshold_percentile", 0.95)
SUMMARY_STATS_COMBOS = config.get("summary_stats_combos", ["QST,F_within_pop"])
CHROMOSOMES = config.get("chromosomes", ["autosomes", "chrX"])
OUTPUT_DIR = config.get("output_dir", "results/evaluate")

# Calculate number of batches
N_NEUTRAL_BATCHES = (NUM_NEUTRAL + BATCH_SIZE - 1) // BATCH_SIZE
N_ADAPTIVE_BATCHES = (NUM_REPEATS + BATCH_SIZE - 1) // BATCH_SIZE

# Rule: Estimate neutral QST for evaluation (per ve_ratio)
rule evaluate_neutral_qst:
    input:
        fst_file=lambda wildcards: config.get(f"fst_{wildcards.chr}", f"data/example/qst_neutral_{wildcards.chr}.txt")
    output:
        neutral_batch=f"{OUTPUT_DIR}/{{chr}}/{{stats}}/neutral_r{{ratio}}_b{{batch}}.RData"
    params:
        ve_ratio="{ratio}",
        batch="{batch}",
        batch_size=BATCH_SIZE,
        num_sim=NUM_SIM,
        summary_stats=lambda wildcards: wildcards.stats.replace("_", ","),
        scripts_dir=SCRIPTS_DIR
    threads: config.get("threads_per_job", 1)
    log:
        f"{OUTPUT_DIR}/logs/{{chr}}_{{stats}}_neutral_r{{ratio}}_b{{batch}}.log"
    shell:
        """
        # Calculate start and end indices for this batch
        START=$(( ({params.batch} - 1) * {params.batch_size} + 1 ))
        END=$(( {params.batch} * {params.batch_size} ))
        
        # Extract FST values for this batch
        FST_BATCH=$(mktemp)
        sed -n "${{START}},${{END}}p" {input.fst_file} > $FST_BATCH
        
        Rscript {params.scripts_dir}/qst_abc_sim.R \
            batch_neutral \
            $FST_BATCH \
            {params.ve_ratio} \
            {output.neutral_batch} \
            {params.num_sim} \
            {params.summary_stats} \
            > {log} 2>&1
        
        rm -f $FST_BATCH
        """

# Rule: Estimate adaptive QST for evaluation
rule evaluate_adaptive_qst:
    input:
        fst_file=lambda wildcards: config.get(f"fst_{wildcards.chr}", f"data/example/qst_neutral_{wildcards.chr}.txt")
    output:
        adaptive_batch=f"{OUTPUT_DIR}/{{chr}}/{{stats}}/adaptive_q{{qst}}_r{{ratio}}_b{{batch}}.RData"
    params:
        adaptive_qst="{qst}",
        ve_ratio="{ratio}",
        batch="{batch}",
        batch_size=BATCH_SIZE,
        num_sim=NUM_SIM,
        summary_stats=lambda wildcards: wildcards.stats.replace("_", ","),
        scripts_dir=SCRIPTS_DIR
    threads: config.get("threads_per_job", 1)
    log:
        f"{OUTPUT_DIR}/logs/{{chr}}_{{stats}}_adaptive_q{{qst}}_r{{ratio}}_b{{batch}}.log"
    shell:
        """
        Rscript {params.scripts_dir}/qst_abc_sim.R \
            batch_evaluate \
            {params.adaptive_qst},{params.ve_ratio},{params.batch_size} \
            {params.ve_ratio} \
            {output.adaptive_batch} \
            {params.num_sim} \
            {params.summary_stats} \
            > {log} 2>&1
        """

# Rule: Aggregate TPR/FPR for one chromosome and summary stats combo
rule aggregate_perf_eval:
    input:
        neutral_batches=lambda wildcards: expand(
            f"{OUTPUT_DIR}/{{chr}}/{{stats}}/neutral_r{{ratio}}_b{{batch}}.RData",
            chr=wildcards.chr,
            stats=wildcards.stats,
            ratio=[str(r) for r in VE_RATIOS],
            batch=range(1, N_NEUTRAL_BATCHES + 1)
        ),
        adaptive_batches=lambda wildcards: expand(
            f"{OUTPUT_DIR}/{{chr}}/{{stats}}/adaptive_q{{qst}}_r{{ratio}}_b{{batch}}.RData",
            chr=wildcards.chr,
            stats=wildcards.stats,
            qst=[str(q) for q in ADAPTIVE_QST_LEVELS],
            ratio=[str(r) for r in VE_RATIOS],
            batch=range(1, N_ADAPTIVE_BATCHES + 1)
        )
    output:
        matrix=f"{OUTPUT_DIR}/{{chr}}/tpr_fpr_matrix_{{chr}}_{{stats}}.csv",
        heatmap=f"{OUTPUT_DIR}/{{chr}}/tpr_fpr_matrix_{{chr}}_{{stats}}_heatmap.pdf"
    params:
        input_dir=f"{OUTPUT_DIR}/{{chr}}/{{stats}}",
        chr="{chr}",
        summary_stats=lambda wildcards: wildcards.stats.replace("_", ","),
        scripts_dir=SCRIPTS_DIR
    log:
        f"{OUTPUT_DIR}/logs/{{chr}}_{{stats}}_aggregate.log"
    shell:
        """
        Rscript {params.scripts_dir}/aggregate_perf_eval.R \
            {params.input_dir} \
            {params.chr} \
            {params.summary_stats} \
            > {log} 2>&1
        """

# Rule: Generate combined model ranking
rule combined_model_ranking:
    input:
        matrices=lambda wildcards: expand(
            f"{OUTPUT_DIR}/{{chr}}/tpr_fpr_matrix_{{chr}}_{{stats}}.csv",
            chr=CHROMOSOMES,
            stats=[s.replace(",", "_") for s in SUMMARY_STATS_COMBOS]
        )
    output:
        ranking=f"{OUTPUT_DIR}/combined_model_ranking.csv"
    params:
        input_dir=OUTPUT_DIR,
        scripts_dir=SCRIPTS_DIR
    log:
        f"{OUTPUT_DIR}/logs/combined_ranking.log"
    shell:
        """
        Rscript {params.scripts_dir}/generate_combined_ranking.R \
            {params.input_dir} \
            {output.ranking} \
            > {log} 2>&1
        """

# Rule: Generate TPR performance plots
rule plot_tpr_performance:
    input:
        ranking=f"{OUTPUT_DIR}/combined_model_ranking.csv"
    output:
        combined_plot=f"{OUTPUT_DIR}/plots/TPR_plot_combined.pdf"
    params:
        input_dir=OUTPUT_DIR,
        output_dir=f"{OUTPUT_DIR}/plots",
        scripts_dir=SCRIPTS_DIR
    log:
        f"{OUTPUT_DIR}/logs/plot_tpr.log"
    shell:
        """
        mkdir -p {params.output_dir}
        Rscript {params.scripts_dir}/plot_tpr_performance.R \
            {params.input_dir} \
            {params.output_dir} \
            > {log} 2>&1
        """
