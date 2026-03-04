# REDQuanTA Module 1: Adaptive QST Detection Rules
#
# Workflow:
#   1. Extract trait IDs from trait_values.csv
#   2. For each trait: prepare_obs_stats -> estimate_trait_qst + estimate_neutral_qst -> aggregate_trait
#   3. Combine all trait results -> qst_results.csv

import pandas as pd
from pathlib import Path

# Get config values
SAMPLE_STRUCTURE = config.get("sample_structure", "data/example/sample_structure.csv")
TRAIT_VALUES = config.get("trait_values", "data/example/trait_values.csv")
NUM_NEUTRAL = config.get("num_neutral", 100)
NUM_SIM = config.get("num_sim", 100000)
BATCH_SIZE = config.get("batch_size", 50)
THRESHOLD_PERCENTILE = config.get("threshold_percentile", 0.95)
SUMMARY_STATS = config.get("summary_stats", "QST,F_within_pop")
CHROMOSOMES = config.get("chromosomes", ["autosomes", "chrX"])
OUTPUT_DIR = config.get("output_dir", "results/detect")

# Calculate number of batches
N_BATCHES = (NUM_NEUTRAL + BATCH_SIZE - 1) // BATCH_SIZE

# Function to get trait IDs from trait_values.csv
def get_trait_ids():
    if os.path.exists(TRAIT_VALUES):
        df = pd.read_csv(TRAIT_VALUES)
        if "trait_id" in df.columns:
            return df["trait_id"].unique().tolist()
    return []

# Get trait IDs at workflow start
TRAIT_IDS = get_trait_ids()

# Rule: Prepare observed statistics for a trait
rule prepare_obs_stats:
    input:
        sample_structure=SAMPLE_STRUCTURE,
        trait_values=TRAIT_VALUES
    output:
        obs_stats=f"{OUTPUT_DIR}/{{trait}}/{{trait}}_obs_stats.RData",
        ext_sd=f"{OUTPUT_DIR}/{{trait}}/{{trait}}_ext_sd.txt"
    params:
        trait_id="{trait}",
        scripts_dir=SCRIPTS_DIR
    log:
        f"{OUTPUT_DIR}/logs/{{trait}}_prepare_obs_stats.log"
    shell:
        """
        Rscript {params.scripts_dir}/prepare_obs_stats.R \
            {input.sample_structure} \
            {input.trait_values} \
            {params.trait_id} \
            {output.obs_stats} \
            {output.ext_sd} \
            > {log} 2>&1
        """

# Rule: Estimate trait QST using ABC
rule estimate_trait_qst:
    input:
        obs_stats=f"{OUTPUT_DIR}/{{trait}}/{{trait}}_obs_stats.RData"
    output:
        qst=f"{OUTPUT_DIR}/{{trait}}/{{trait}}_trait_qst.RData"
    params:
        num_sim=NUM_SIM,
        summary_stats=SUMMARY_STATS,
        scripts_dir=SCRIPTS_DIR
    threads: config.get("threads_per_job", 1)
    log:
        f"{OUTPUT_DIR}/logs/{{trait}}_estimate_trait_qst.log"
    shell:
        """
        Rscript {params.scripts_dir}/qst_abc_sim.R \
            trait \
            {input.obs_stats} \
            NA \
            {output.qst} \
            {params.num_sim} \
            {params.summary_stats} \
            > {log} 2>&1
        """

# Rule: Estimate neutral QST for a batch of FST values
rule estimate_neutral_qst:
    input:
        ext_sd=f"{OUTPUT_DIR}/{{trait}}/{{trait}}_ext_sd.txt",
        fst_file=lambda wildcards: config.get(f"fst_{wildcards.chr}", f"data/example/qst_neutral_{wildcards.chr}.txt")
    output:
        neutral_batch=f"{OUTPUT_DIR}/{{trait}}/{{chr}}/neutral_batch_{{batch}}.RData"
    params:
        batch="{batch}",
        batch_size=BATCH_SIZE,
        num_sim=NUM_SIM,
        summary_stats=SUMMARY_STATS,
        scripts_dir=SCRIPTS_DIR
    threads: config.get("threads_per_job", 1)
    log:
        f"{OUTPUT_DIR}/logs/{{trait}}_{{chr}}_neutral_batch_{{batch}}.log"
    shell:
        """
        # Calculate start and end indices for this batch
        START=$(( ({params.batch} - 1) * {params.batch_size} + 1 ))
        END=$(( {params.batch} * {params.batch_size} ))
        
        # Extract FST values for this batch
        FST_BATCH=$(mktemp)
        sed -n "${{START}},${{END}}p" {input.fst_file} > $FST_BATCH
        
        # Read ext_sd
        EXT_SD=$(cat {input.ext_sd})
        
        Rscript {params.scripts_dir}/qst_abc_sim.R \
            batch_neutral \
            $FST_BATCH \
            $EXT_SD \
            {output.neutral_batch} \
            {params.num_sim} \
            {params.summary_stats} \
            > {log} 2>&1
        
        rm -f $FST_BATCH
        """

# Rule: Aggregate results for a single trait
rule aggregate_trait:
    input:
        trait_qst=f"{OUTPUT_DIR}/{{trait}}/{{trait}}_trait_qst.RData",
        neutral_batches=lambda wildcards: expand(
            f"{OUTPUT_DIR}/{{trait}}/{{chr}}/neutral_batch_{{batch}}.RData",
            trait=wildcards.trait,
            chr=CHROMOSOMES,
            batch=range(1, N_BATCHES + 1)
        )
    output:
        result=f"{OUTPUT_DIR}/{{trait}}/{{trait}}_result.csv"
    params:
        trait_id="{trait}",
        threshold_percentile=THRESHOLD_PERCENTILE,
        scripts_dir=SCRIPTS_DIR
    log:
        f"{OUTPUT_DIR}/logs/{{trait}}_aggregate.log"
    shell:
        """
        Rscript {params.scripts_dir}/aggregate_qst.R \
            {params.trait_id} \
            {input.trait_qst} \
            {OUTPUT_DIR}/{params.trait_id} \
            {output.result} \
            {params.threshold_percentile} \
            > {log} 2>&1
        """

# Rule: Combine all trait results
rule aggregate_all_traits:
    input:
        results=lambda wildcards: expand(
            f"{OUTPUT_DIR}/{{trait}}/{{trait}}_result.csv",
            trait=TRAIT_IDS
        ) if TRAIT_IDS else []
    output:
        combined=f"{OUTPUT_DIR}/qst_results.csv"
    log:
        f"{OUTPUT_DIR}/logs/aggregate_all.log"
    run:
        import pandas as pd
        
        if not input.results:
            # Create empty result file if no traits
            pd.DataFrame(columns=["trait_id", "chr", "QST", "threshold_percentile", 
                                  "threshold_value", "adaptive"]).to_csv(output.combined, index=False)
        else:
            dfs = [pd.read_csv(f) for f in input.results]
            combined = pd.concat(dfs, ignore_index=True)
            combined.to_csv(output.combined, index=False)
        
        print(f"Combined {len(input.results)} trait results into {output.combined}")
