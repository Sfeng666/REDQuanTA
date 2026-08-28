#!/usr/bin/env python3
"""Generate HTCondor DAGs for fast all-combo QST performance evaluation (REDQuanTEA).

Uses workflow/scripts/{qst_abc_sim.R, abc_fast_estimators.R,
abc_boundary_helpers.R, perf_eval_all_combos_fast_cell_v2.R} and
aggregate_perf_eval_multicombo_fast.R on the submit node after completion.
"""

import argparse
import os
import shutil
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
HTCONDOR_DIR = SCRIPT_DIR.parent
BASE_DIR = HTCONDOR_DIR.parent
WORKFLOW_SCRIPTS_DIR = BASE_DIR / "workflow" / "scripts"
INPUT_DIR_DEFAULT = BASE_DIR / "data" / "example"
R_ENV_DEFAULT = HTCONDOR_DIR / "env" / "r_env.tar.gz"

RSCRIPT = os.environ.get("RSCRIPT", shutil.which("Rscript") or "Rscript")

DEFAULT_ADAPTIVE_QST = [round(0.50 + i * 0.05, 2) for i in range(11)]
DEFAULT_VE_RATIOS = [0.01, 0.1, 1.0, 10.0, 100.0]


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


def write_fst_batches(values, output_dir, batch_size):
    batches = []
    for i in range(0, len(values), batch_size):
        batch = values[i : i + batch_size]
        path = output_dir / f"fst_batch_{len(batches) + 1}.txt"
        with open(path, "w", encoding="utf-8") as f:
            for value in batch:
                f.write(f"{value}\n")
        batches.append((path, len(batch), i + 1))
    return batches


def count_combos(path):
    with open(path, "r", encoding="utf-8") as f:
        return sum(1 for line in f if line.strip())


def aggregate_script(output_dir, chr_type, adaptive_qst, ve_ratios, threshold, combos_file):
    script = output_dir / f"aggregate_perf_eval_{chr_type}.sh"
    result_file = output_dir / f"tpr_fpr_matrix_{chr_type}.csv"
    agg_r = WORKFLOW_SCRIPTS_DIR / "aggregate_perf_eval_multicombo_fast.R"
    with open(script, "w", encoding="utf-8") as f:
        f.write("#!/usr/bin/env bash\n")
        f.write("set -euo pipefail\n")
        f.write("export PERF_EVAL_DETECTOR_STATS=qst,p_qst_gt_fst95\n")
        f.write(f'"{RSCRIPT}" "{agg_r}" \\\n')
        f.write(f'  "{output_dir}" \\\n')
        f.write(f'  "{chr_type}" \\\n')
        f.write(f'  "{",".join(map(str, adaptive_qst))}" \\\n')
        f.write(f'  "{",".join(map(str, ve_ratios))}" \\\n')
        f.write(f'  "{threshold}" \\\n')
        f.write(f'  "{result_file}" \\\n')
        f.write(f'  "{combos_file}"\n')
    os.chmod(script, 0o755)
    return script


def generate_dag(chr_type, args, output_root):
    output_dir = (output_root / chr_type).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    combos_file = Path(args.summary_stats).resolve()
    combos_base = combos_file.name
    code_dir = Path(args.code_dir).resolve() if args.code_dir else WORKFLOW_SCRIPTS_DIR.resolve()

    fst_dir = Path(args.fst_input_dir)
    fst_file = fst_dir / ("qst_neutral_chrX.txt" if chr_type == "chrX" else "qst_neutral_autosomes.txt")
    fst_values = read_fst_values(fst_file, args.num_neutral)
    batch_size = min(args.batch_size, args.num_neutral, args.num_repeats)
    adaptive_qst = parse_float_list(args.adaptive_qst, DEFAULT_ADAPTIVE_QST)
    ve_ratios = parse_float_list(args.ve_ratios, DEFAULT_VE_RATIOS)
    env_vars = (
        f"SIM_NUM_POP={args.num_pop} SIM_NUM_IND={args.num_ind} SIM_NUM_REP={args.num_rep} "
        f"QST_PERF_EVAL_CHR={chr_type} "
        f"QST_ABC_METHOD=loclinear QST_ABC_TOL_NUMERATOR={args.abc_tol_numerator} "
        f"FAST_ABC_ESTIMATOR={args.estimator} FAST_ABC_COMBO_CHUNK_SIZE={args.combo_chunk_size} "
        f"POC_P2_MC_CORES={args.mc_cores} QST_ABC_BATCH_SIZE=50000 "
        f"FLOOR_POLICY={args.floor_policy} FLOOR_ALPHA={args.floor_alpha}"
    )
    sub_path = SCRIPT_DIR / args.sub_name
    if "://" in args.r_env_tarball:
        r_env = args.r_env_tarball
    else:
        r_env = str(Path(args.r_env_tarball).resolve())

    dag_path = output_root / f"perf_eval_{chr_type}.dag"
    all_jobs = []
    with open(dag_path, "w", encoding="utf-8") as f:
        f.write("# HTCondor DAG for Design Module QST performance evaluation (REDQuanTEA)\n")
        f.write(f"# Chromosome: {chr_type}\n")
        f.write(f"# Combos: {count_combos(combos_file)} from {combos_file}\n")
        f.write(f"# Submit template: {sub_path}\n")
        f.write(f"CONFIG {SCRIPT_DIR / 'unlimited.config'}\n\n")

        for ratio_idx, ve_ratio in enumerate(ve_ratios, 1):
            ratio_dir = output_dir / f"neutral_ratio_{ratio_idx}"
            ratio_dir.mkdir(parents=True, exist_ok=True)
            for batch_idx, (batch_file, _n, start_idx) in enumerate(
                write_fst_batches(fst_values, ratio_dir, batch_size), 1
            ):
                job = f"neutral_r{ratio_idx}_b{batch_idx}"
                out = f"neutral_batch_{batch_idx}.RData"
                worker_args = f'{batch_file.name} {start_idx} {ve_ratio} {combos_base} {out} {args.num_sim}'
                f.write(f"JOB {job} {sub_path}\n")
                f.write(
                    f'VARS {job} mode="neutral" worker_args="{worker_args}" output_file="{out}" '
                    f'outdir="{ratio_dir}" dir_scripts="{SCRIPT_DIR}" dir_code="{code_dir}" '
                    f'combos_file="{combos_file}" r_env_tarball="{r_env}" '
                    f'extra_input_files=", {batch_file.resolve()}" request_memory="{args.request_memory}" job="{job}"\n'
                )
                f.write(f'VARS {job} environment="{env_vars}"\n')
                all_jobs.append(job)
            f.write("\n")

        for qst in adaptive_qst:
            qst_str = f"{qst:.2f}".replace(".", "_")
            for ratio_idx, ve_ratio in enumerate(ve_ratios, 1):
                cond_dir = output_dir / f"adaptive_q{qst_str}_r{ratio_idx}"
                cond_dir.mkdir(parents=True, exist_ok=True)
                n_batches = (args.num_repeats + batch_size - 1) // batch_size
                for batch_idx in range(n_batches):
                    start_id = batch_idx * batch_size + 1
                    n_in_batch = min(batch_size, args.num_repeats - batch_idx * batch_size)
                    job = f"adaptive_q{qst_str}_r{ratio_idx}_b{batch_idx + 1}"
                    out = f"adaptive_batch_{batch_idx + 1}.RData"
                    worker_args = f"{start_id} {n_in_batch} {qst} {ve_ratio} {combos_base} {out} {args.num_sim}"
                    f.write(f"JOB {job} {sub_path}\n")
                    f.write(
                        f'VARS {job} mode="adaptive" worker_args="{worker_args}" output_file="{out}" '
                        f'outdir="{cond_dir}" dir_scripts="{SCRIPT_DIR}" dir_code="{code_dir}" '
                        f'combos_file="{combos_file}" r_env_tarball="{r_env}" '
                        f'extra_input_files="" request_memory="{args.request_memory}" job="{job}"\n'
                    )
                    f.write(f'VARS {job} environment="{env_vars}"\n')
                    all_jobs.append(job)
                f.write("\n")

        noop = f"noop_perf_eval_{chr_type}"
        f.write(f"JOB {noop} {SCRIPT_DIR / 'noop.sub'} NOOP\n")
        f.write(f"PARENT {' '.join(all_jobs)} CHILD {noop}\n")
        f.write(f"SCRIPT POST {noop} {aggregate_script(output_dir, chr_type, adaptive_qst, ve_ratios, args.threshold_percentile, combos_file.resolve())}\n")

    return dag_path, len(all_jobs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-stats", required=True, help="Tab-separated combo file")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--chr", choices=("autosomes", "chrX", "both"), default="both")
    parser.add_argument("--adaptive-qst")
    parser.add_argument("--ve-ratios")
    parser.add_argument("--num-repeats", type=int, default=10000)
    parser.add_argument("--num-neutral", type=int, default=10000)
    parser.add_argument("--num-sim", type=int, default=100000)
    parser.add_argument("--batch-size", type=int, default=1000)
    parser.add_argument("--threshold-percentile", type=float, default=0.95)
    parser.add_argument("--abc-tol-numerator", type=int, default=50)
    parser.add_argument("--num-pop", type=int, default=2)
    parser.add_argument("--num-ind", type=int, default=6)
    parser.add_argument("--num-rep", type=int, default=3)
    parser.add_argument("--estimator", choices=("loclinear", "rejection"), default="loclinear")
    parser.add_argument("--combo-chunk-size", type=int, default=64)
    parser.add_argument("--request-memory", default="4GB")
    parser.add_argument("--mc-cores", type=int, default=4, help="POC_P2_MC_CORES for worker (use 1 on OSPool)")
    parser.add_argument("--code-dir", help="Directory with qst_abc_sim.R (default: workflow/scripts)")
    parser.add_argument("--fst-input-dir", default=str(INPUT_DIR_DEFAULT), help="Directory with qst_neutral_*.txt")
    parser.add_argument("--r-env-tarball", default=str(R_ENV_DEFAULT.resolve()))
    parser.add_argument("--floor-policy", default="ridge_floor",
                        help="Variance floor (default: ridge_floor)")
    parser.add_argument("--floor-alpha", default="0.1",
                        help="Ridge-floor alpha (default: 0.1)")
    parser.add_argument(
        "--sub-name",
        default="perf_eval_fast.sub",
        help="Submit template under htcondor/scripts/",
    )
    args = parser.parse_args()

    output_root = Path(args.output_dir).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    chroms = ["autosomes", "chrX"] if args.chr == "both" else [args.chr]
    for chr_type in chroms:
        dag_path, n_jobs = generate_dag(chr_type, args, output_root)
        print(f"Wrote {dag_path} with {n_jobs} jobs")


if __name__ == "__main__":
    main()
