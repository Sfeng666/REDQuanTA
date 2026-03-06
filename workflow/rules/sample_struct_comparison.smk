# REDQuanTA: Sample Structure Comparison Rules
SAMPLE_STRUCTS = config.get("sample_structures", [])
COMPARISON_STATS = config.get("comparison_stats",
    config.get("summary_stats_combos", ["QST,F_within_pop"])[0])
COMPARISON_STATS_SAFE = COMPARISON_STATS.replace(",", "_")
SS_BASE = f"{OUTPUT_DIR}/sample_struct_comparison_{COMPARISON_STATS_SAFE}"

STRUCT_MAP = {}
STRUCT_IDS = []
for _s in SAMPLE_STRUCTS:
    _sid = f"n2_i{_s['n_ind']}_r{_s['n_rep']}"
    STRUCT_IDS.append(_sid)
    STRUCT_MAP[_sid] = _s

if SAMPLE_STRUCTS:
    rule ss_neutral_qst:
        input:
            fst_file=lambda wc: config.get(f"fst_{wc.chr}", f"data/example/qst_neutral_{wc.chr}.txt")
        output:
            batch=f"{SS_BASE}/{{struct}}/{{chr}}/neutral_r{{ratio}}_b{{batch}}.RData"
        params:
            ve_ratio="{ratio}", batch="{batch}",
            batch_size=BATCH_SIZE, num_sim=NUM_SIM,
            summary_stats=COMPARISON_STATS,
            n_ind=lambda wc: STRUCT_MAP[wc.struct]['n_ind'],
            n_rep=lambda wc: STRUCT_MAP[wc.struct]['n_rep'],
            scripts_dir=SCRIPTS_DIR
        threads: config.get("threads_per_job", 1)
        log: f"{OUTPUT_DIR}/logs/ss_{{struct}}_{{chr}}_neutral_r{{ratio}}_b{{batch}}.log"
        shell:
            """
            mkdir -p $(dirname {output.batch})
            START=$(( ({params.batch} - 1) * {params.batch_size} + 1 ))
            END=$(( {params.batch} * {params.batch_size} ))
            FST_BATCH=$(mktemp)
            sed -n "${{START}},${{END}}p" {input.fst_file} > $FST_BATCH
            SIM_NUM_IND={params.n_ind} SIM_NUM_REP={params.n_rep} \
            Rscript {params.scripts_dir}/qst_abc_sim.R \
                batch_neutral $FST_BATCH {params.ve_ratio} \
                {output.batch} {params.num_sim} {params.summary_stats} \
                > {log} 2>&1
            rm -f $FST_BATCH
            """

    rule ss_adaptive_qst:
        input:
            fst_file=lambda wc: config.get(f"fst_{wc.chr}", f"data/example/qst_neutral_{wc.chr}.txt")
        output:
            batch=f"{SS_BASE}/{{struct}}/{{chr}}/adaptive_q{{qst}}_r{{ratio}}_b{{batch}}.RData"
        params:
            adaptive_qst="{qst}", ve_ratio="{ratio}",
            batch="{batch}", batch_size=BATCH_SIZE,
            num_sim=NUM_SIM, summary_stats=COMPARISON_STATS,
            n_ind=lambda wc: STRUCT_MAP[wc.struct]['n_ind'],
            n_rep=lambda wc: STRUCT_MAP[wc.struct]['n_rep'],
            scripts_dir=SCRIPTS_DIR
        threads: config.get("threads_per_job", 1)
        log: f"{OUTPUT_DIR}/logs/ss_{{struct}}_{{chr}}_adaptive_q{{qst}}_r{{ratio}}_b{{batch}}.log"
        shell:
            """
            mkdir -p $(dirname {output.batch})
            START_ID=$(( ({params.batch} - 1) * {params.batch_size} + 1 ))
            SIM_NUM_IND={params.n_ind} SIM_NUM_REP={params.n_rep} \
            Rscript {params.scripts_dir}/qst_abc_sim.R \
                batch_evaluate \
                "${{START_ID}}_{params.batch_size}_{params.adaptive_qst}" \
                {params.ve_ratio} {output.batch} \
                {params.num_sim} {params.summary_stats} \
                > {log} 2>&1
            """

    rule ss_aggregate:
        input:
            neutral=lambda wc: expand(
                f"{SS_BASE}/{wc.struct}/{wc.chr}/neutral_r{{ratio}}_b{{batch}}.RData",
                ratio=[str(r) for r in VE_RATIOS],
                batch=range(1, N_NEUTRAL_BATCHES + 1)),
            adaptive=lambda wc: expand(
                f"{SS_BASE}/{wc.struct}/{wc.chr}/adaptive_q{{qst}}_r{{ratio}}_b{{batch}}.RData",
                qst=[str(q) for q in ADAPTIVE_QST_LEVELS],
                ratio=[str(r) for r in VE_RATIOS],
                batch=range(1, N_ADAPTIVE_BATCHES + 1))
        output:
            matrix=f"{SS_BASE}/{{struct}}/{{chr}}/tpr_fpr_matrix_{{chr}}.csv"
        params:
            input_dir=f"{SS_BASE}/{{struct}}/{{chr}}",
            chr="{chr}",
            adaptive_qst=",".join([str(q) for q in ADAPTIVE_QST_LEVELS]),
            ve_ratios=",".join([str(r) for r in VE_RATIOS]),
            threshold_percentile=THRESHOLD_PERCENTILE,
            scripts_dir=SCRIPTS_DIR
        log: f"{OUTPUT_DIR}/logs/ss_{{struct}}_{{chr}}_aggregate.log"
        shell:
            """
            Rscript {params.scripts_dir}/aggregate_perf_eval.R \
                {params.input_dir} {params.chr} \
                {params.adaptive_qst} {params.ve_ratios} \
                {params.threshold_percentile} {output.matrix} \
                > {log} 2>&1
            """

    rule ss_plot:
        input:
            matrices=expand(
                f"{SS_BASE}/{{struct}}/{{chr}}/tpr_fpr_matrix_{{chr}}.csv",
                struct=STRUCT_IDS, chr=CHROMOSOMES)
        output:
            done=f"{SS_BASE}/plots/done.txt"
        params:
            base_dir=SS_BASE,
            output_dir=f"{SS_BASE}/plots",
            summary_stats=COMPARISON_STATS,
            scripts_dir=SCRIPTS_DIR
        log: f"{OUTPUT_DIR}/logs/ss_plot.log"
        shell:
            """
            mkdir -p {params.output_dir}
            Rscript {params.scripts_dir}/plot_sample_structure_comparison.R \
                {params.base_dir} {params.output_dir} \
                {params.summary_stats} \
                > {log} 2>&1
            touch {output.done}
            """
