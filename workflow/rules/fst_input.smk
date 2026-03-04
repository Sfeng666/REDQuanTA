# REDQuanTA FST Input Rules
#
# Three modes for obtaining neutral FST distribution:
#   1. direct: User provides FST file paths in config (default)
#   2. from_vcf: Generate FST from VCF file
#   3. from_simulation: Generate FST from ms coalescent simulation

import os
from pathlib import Path

FST_INPUT_MODE = config.get("fst_input_mode", "direct")
FST_OUTPUT_DIR = config.get("output_dir", "results") + "/fst_input"

# Rule: Generate FST from VCF (mode: from_vcf)
rule fst_from_vcf:
    input:
        vcf=config.get("vcf_file", "")
    output:
        fst_autosomes=f"{FST_OUTPUT_DIR}/qst_neutral_autosomes.txt",
        fst_chrX=f"{FST_OUTPUT_DIR}/qst_neutral_chrX.txt"
    params:
        populations=config.get("vcf_populations", ""),
        scripts_dir=SCRIPTS_DIR
    log:
        f"{FST_OUTPUT_DIR}/fst_from_vcf.log"
    shell:
        """
        python {params.scripts_dir}/fst_from_vcf.py \
            --vcf {input.vcf} \
            --populations {params.populations} \
            --output-autosomes {output.fst_autosomes} \
            --output-chrX {output.fst_chrX} \
            > {log} 2>&1
        """

# Rule: Generate FST from ms simulation (mode: from_simulation)
rule fst_from_simulation:
    input:
        ms_config=config.get("ms_config", "")
    output:
        fst_autosomes=f"{FST_OUTPUT_DIR}/qst_neutral_autosomes.txt",
        fst_chrX=f"{FST_OUTPUT_DIR}/qst_neutral_chrX.txt"
    params:
        scripts_dir=SCRIPTS_DIR
    log:
        f"{FST_OUTPUT_DIR}/fst_from_simulation.log"
    shell:
        """
        python {params.scripts_dir}/fst_from_simulation.py \
            --config {input.ms_config} \
            --output-autosomes {output.fst_autosomes} \
            --output-chrX {output.fst_chrX} \
            > {log} 2>&1
        """

# Helper function to determine FST file path based on mode
def get_fst_file(chr_type):
    mode = config.get("fst_input_mode", "direct")
    
    if mode == "direct":
        return config.get(f"fst_{chr_type}", f"data/example/qst_neutral_{chr_type}.txt")
    elif mode in ["from_vcf", "from_simulation"]:
        return f"{FST_OUTPUT_DIR}/qst_neutral_{chr_type}.txt"
    else:
        raise ValueError(f"Unknown fst_input_mode: {mode}")

# Rule: Prepare FST input (dispatcher based on mode)
rule prepare_fst_input:
    input:
        fst_autosomes=lambda wildcards: get_fst_file("autosomes"),
        fst_chrX=lambda wildcards: get_fst_file("chrX")
    output:
        touch(f"{FST_OUTPUT_DIR}/.fst_ready")
    shell:
        """
        echo "FST input ready:"
        echo "  Autosomes: {input.fst_autosomes}"
        echo "  chrX: {input.fst_chrX}"
        """
