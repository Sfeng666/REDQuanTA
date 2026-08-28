#!/usr/bin/env python3
"""Generate HTCondor DAGs for aligned Design Module ANOVA perf-eval (mega-batch jobs)."""

import argparse
import json
import os
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
HTCONDOR_DIR = SCRIPT_DIR.parent
BASE_DIR = HTCONDOR_DIR.parent
WORKFLOW_SCRIPTS_DIR = BASE_DIR / "workflow" / "scripts"
INPUT_DIR_DEFAULT = BASE_DIR / "data" / "example"
R_ENV_DEFAULT = HTCONDOR_DIR / "env" / "r_env.tar.gz"

RSCRIPT = os.environ.get("RSCRIPT", "Rscript")

DEFAULT_ADAPTIVE_QST = [round(0.50 + i * 0.05, 2) for i in range(11)]
DEFAULT_VE_RATIOS = [0.01, 0.1, 1.0, 10.0, 100.0]


def neutral_output_dirs(n_ratios=None):
    n_ratios = n_ratios or len(DEFAULT_VE_RATIOS)
    return ",".join(f"neutral_ratio_{i}" for i in range(1, n_ratios + 1))


def adaptive_output_dirs(adaptive_qst, n_ratios=None):
    n_ratios = n_ratios or len(DEFAULT_VE_RATIOS)
    dirs = []
    for qst in adaptive_qst:
        qst_str = f"{qst:.2f}".replace(".", "_")
        for ri in range(1, n_ratios + 1):
            dirs.append(f"adaptive_q{qst_str}_r{ri}")
    return ",".join(dirs)


def adaptive_q_output_dirs(qst, n_ratios=None):
    n_ratios = n_ratios or len(DEFAULT_VE_RATIOS)
    qst_str = f"{qst:.2f}".replace(".", "_")
    return ",".join(f"adaptive_q{qst_str}_r{ri}" for ri in range(1, n_ratios + 1))


def parse_float_list(value, default):
    if value is None:
        return default
    return [float(x) for x in value.split(",") if x.strip()]


def read_fst_values(path, n):
    values = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if len(values) >= n:
                break
            line = line.strip()
            if line:
                values.append(float(line))
    return values


def write_combo_file(output_dir, combo_name):
    path = output_dir / f"combo_{combo_name}.txt"
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{combo_name}\n")
    return path


def aggregate_script(output_dir, chr_type, adaptive_qst, ve_ratios, threshold, combo_file):
    script = output_dir / f"aggregate_perf_eval_{chr_type}.sh"
    result_file = output_dir / f"tpr_fpr_matrix_{chr_type}.csv"
    agg_r = WORKFLOW_SCRIPTS_DIR / "aggregate_perf_eval_multicombo_fast.R"
    with open(script, "w", encoding="utf-8") as f:
        f.write("#!/usr/bin/env bash\n")
        f.write("set -euo pipefail\n")
        f.write(f'"{RSCRIPT}" "{agg_r}" \\\n')
        f.write(f'  "{output_dir}" \\\n')
        f.write(f'  "{chr_type}" \\\n')
        f.write(f'  "{",".join(map(str, adaptive_qst))}" \\\n')
        f.write(f'  "{",".join(map(str, ve_ratios))}" \\\n')
        f.write(f'  "{threshold}" \\\n')
        f.write(f'  "{result_file}" \\\n')
        f.write(f'  "{combo_file}"\n')
    os.chmod(script, 0o755)
    return script


def adaptive_split_strategy(args, num_trials, n_qst, n_ve):
    """Return 'mega' or 'by_qst' based on benchmark seconds/trial."""
    sec_per_trial = args.sec_per_trial
    if sec_per_trial is None and args.benchmark_json:
        p = Path(args.benchmark_json)
        if p.exists():
            data = json.loads(p.read_text(encoding="utf-8"))
            sec_per_trial = data.get("sec_per_trial_adaptive")
    if sec_per_trial is None:
        sec_per_trial = float(os.environ.get("ANOVA_SEC_PER_TRIAL", "0"))
    if sec_per_trial <= 0:
        return "mega"
    total_trials = n_qst * n_ve * num_trials
    est_minutes = (total_trials * sec_per_trial) / max(args.request_cpus, 1) / 60.0
    return "by_qst" if est_minutes > args.max_wall_minutes else "mega"


def generate_dag(chr_type, args, output_root):
    output_dir = (output_root / chr_type).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    code_dir = Path(args.code_dir).resolve() if args.code_dir else WORKFLOW_SCRIPTS_DIR.resolve()

    combo_name = "ANOVA_norep" if args.method == "anova_norep" else "ANOVA"
    combo_file = write_combo_file(output_dir, combo_name)

    fst_dir = Path(args.fst_input_dir)
    fst_file = fst_dir / ("qst_neutral_chrX.txt" if chr_type == "chrX" else "qst_neutral_autosomes.txt")
    fst_values = read_fst_values(fst_file, args.num_trials)
    fst_local = output_dir / "fst_all.txt"
    with open(fst_local, "w", encoding="utf-8") as f:
        for v in fst_values:
            f.write(f"{v}\n")

    adaptive_qst = parse_float_list(args.adaptive_qst, DEFAULT_ADAPTIVE_QST)
    ve_ratios = parse_float_list(args.ve_ratios, DEFAULT_VE_RATIOS)
    split = adaptive_split_strategy(args, args.num_trials, len(adaptive_qst), len(ve_ratios))

    norep_flag = "ANOVA_NOREP=1" if args.method == "anova_norep" else "ANOVA_NOREP=0"
    env_vars = (
        f"SIM_NUM_POP={args.num_pop} SIM_NUM_IND={args.num_ind} SIM_NUM_REP={args.num_rep} "
        f"POC_P2_MC_CORES={args.request_cpus} {norep_flag}"
    )

    sub_path = SCRIPT_DIR / args.sub_name
    if "://" in args.r_env_tarball:
        r_env = args.r_env_tarball
    else:
        r_env = str(Path(args.r_env_tarball).resolve())

    dag_path = output_root / f"perf_eval_anova_{chr_type}.dag"
    all_jobs = []
    priority = args.priority

    with open(dag_path, "w", encoding="utf-8") as f:
        f.write("# HTCondor DAG for aligned Design Module ANOVA perf-eval (mega-batch)\n")
        f.write(f"# Chromosome: {chr_type} method: {args.method} adaptive_split: {split}\n")
        f.write(f"CONFIG {SCRIPT_DIR / 'unlimited.config'}\n\n")

        job = f"neutral_{chr_type}"
        worker_args = f'neutral_chr fst_all.txt . {args.num_trials} {args.batch_size}'
        out_files = neutral_output_dirs(len(ve_ratios))
        f.write(f"JOB {job} {sub_path}\n")
        f.write(
            f'VARS {job} worker_args="{worker_args}" output_files="{out_files}" output_remaps="" '
            f'outdir="{output_dir}" dir_scripts="{SCRIPT_DIR}" dir_code="{code_dir}" '
            f'r_env_tarball="{r_env}" extra_input_files=", {fst_local.resolve()}" '
            f'request_memory="{args.request_memory}" request_cpus="{args.request_cpus}" job="{job}"\n'
        )
        f.write(f'VARS {job} environment="{env_vars}"\n')
        if priority is not None:
            f.write(f"PRIORITY {job} {priority}\n")
        all_jobs.append(job)
        f.write("\n")

        if split == "mega":
            job = f"adaptive_{chr_type}"
            worker_args = f"adaptive_chr . {args.num_trials} {args.batch_size}"
            out_files = adaptive_output_dirs(adaptive_qst, len(ve_ratios))
            f.write(f"JOB {job} {sub_path}\n")
            f.write(
                f'VARS {job} worker_args="{worker_args}" output_files="{out_files}" output_remaps="" '
                f'outdir="{output_dir}" dir_scripts="{SCRIPT_DIR}" dir_code="{code_dir}" '
                f'r_env_tarball="{r_env}" extra_input_files="" '
                f'request_memory="{args.request_memory}" request_cpus="{args.request_cpus}" job="{job}"\n'
            )
            f.write(f'VARS {job} environment="{env_vars}"\n')
            if priority is not None:
                f.write(f"PRIORITY {job} {priority}\n")
            all_jobs.append(job)
            f.write("\n")
        else:
            for qst in adaptive_qst:
                qst_str = f"{qst:.2f}".replace(".", "_")
                job = f"adaptive_q{qst_str}_{chr_type}"
                worker_args = f"adaptive_q {qst} . {args.num_trials} {args.batch_size}"
                out_files = adaptive_q_output_dirs(qst, len(ve_ratios))
                f.write(f"JOB {job} {sub_path}\n")
                f.write(
                    f'VARS {job} worker_args="{worker_args}" output_files="{out_files}" output_remaps="" '
                    f'outdir="{output_dir}" dir_scripts="{SCRIPT_DIR}" dir_code="{code_dir}" '
                    f'r_env_tarball="{r_env}" extra_input_files="" '
                    f'request_memory="{args.request_memory}" request_cpus="{args.request_cpus}" job="{job}"\n'
                )
                f.write(f'VARS {job} environment="{env_vars}"\n')
                if priority is not None:
                    f.write(f"PRIORITY {job} {priority}\n")
                all_jobs.append(job)
            f.write("\n")

        noop = f"noop_anova_{chr_type}"
        f.write(f"JOB {noop} {SCRIPT_DIR / 'noop.sub'} NOOP\n")
        f.write(f"PARENT {' '.join(all_jobs)} CHILD {noop}\n")
        agg = aggregate_script(
            output_dir, chr_type, adaptive_qst, ve_ratios, args.threshold_percentile, combo_file.resolve()
        )
        f.write(f"SCRIPT POST {noop} {agg}\n")
        if priority is not None:
            f.write(f"PRIORITY {noop} {priority}\n")

    return dag_path, len(all_jobs), split


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", choices=("anova", "anova_norep"), required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--chr", choices=("autosomes", "chrX", "both"), default="both")
    parser.add_argument("--adaptive-qst")
    parser.add_argument("--ve-ratios")
    parser.add_argument("--num-trials", type=int, default=10000)
    parser.add_argument("--batch-size", type=int, default=1000)
    parser.add_argument("--max-wall-minutes", type=float, default=30.0)
    parser.add_argument("--sec-per-trial", type=float, help="Benchmark adaptive sec/trial (with-rep or norep)")
    parser.add_argument("--benchmark-json", help="JSON with sec_per_trial_adaptive from local benchmark")
    parser.add_argument("--threshold-percentile", type=float, default=0.95)
    parser.add_argument("--num-pop", type=int, default=2)
    parser.add_argument("--num-ind", type=int, default=6)
    parser.add_argument("--num-rep", type=int, default=3)
    parser.add_argument("--request-cpus", type=int, default=4)
    parser.add_argument("--request-memory", default="2GB")
    parser.add_argument("--priority", type=int, default=None)
    parser.add_argument("--code-dir", help="Directory with qst_abc_sim.R")
    parser.add_argument("--fst-input-dir", default=str(INPUT_DIR_DEFAULT))
    parser.add_argument("--r-env-tarball", default=str(R_ENV_DEFAULT.resolve()))
    parser.add_argument("--sub-name", default="perf_eval_anova_megabatch_flock.sub")
    args = parser.parse_args()

    output_root = Path(args.output_dir).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    chroms = ["autosomes", "chrX"] if args.chr == "both" else [args.chr]
    for chr_type in chroms:
        dag_path, n_jobs, split = generate_dag(chr_type, args, output_root)
        print(f"Wrote {dag_path} with {n_jobs} compute jobs (adaptive_split={split})")


if __name__ == "__main__":
    main()
