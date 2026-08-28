#!/usr/bin/env python3
"""Generate HTCondor DAGs for sample-structure comparison (single combo, v2 fast ABC).

Wraps prepare_perf_eval_all_combos_fast_dag.py with sample-structure defaults matching
validation_sample_struct_QST_ratioVbetweenVtotal, but uses megabatches (batch_size=2500)
so each job runs ~60–90 min instead of thousands of short jobs.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PREP = SCRIPT_DIR / "prepare_design_dag.py"

DEFAULT_STRUCTURES = [
    (4, 3),
    (6, 3),
    (8, 3),
    (10, 3),
    (12, 3),
    (14, 3),
    (6, 2),
    (9, 2),
    (12, 2),
    (15, 2),
    (18, 2),
    (21, 2),
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", required=True, help="validation_sample_struct_original root")
    parser.add_argument("--combo-file", required=True, help="Single-line tab-separated combo file")
    parser.add_argument("--code-dir", required=True)
    parser.add_argument("--fst-input-dir", required=True)
    parser.add_argument("--r-env-tarball", required=True)
    parser.add_argument("--sub-name", default="perf_eval_fast.sub")
    parser.add_argument("--batch-size", type=int, default=2500)
    parser.add_argument("--num-neutral", type=int, default=10000)
    parser.add_argument("--num-repeats", type=int, default=10000)
    parser.add_argument("--num-sim", type=int, default=100000)
    parser.add_argument("--adaptive-qst", default="0.8")
    parser.add_argument("--ve-ratios", default="0.1,1.0,10.0")
    parser.add_argument("--structures", help="Optional: '4:3,6:3,...' overrides default list")
    parser.add_argument("--chr", choices=("autosomes", "chrX", "both"), default="both")
    args = parser.parse_args()

    output_root = Path(args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    if args.structures:
        structures = []
        for item in args.structures.split(","):
            n_ind, n_rep = item.strip().split(":")
            structures.append((int(n_ind), int(n_rep)))
    else:
        structures = DEFAULT_STRUCTURES

    for n_ind, n_rep in structures:
        struct_id = f"n2_i{n_ind}_r{n_rep}"
        struct_dir = output_root / struct_id
        cmd = [
            sys.executable,
            str(PREP),
            "--summary-stats",
            str(Path(args.combo_file).resolve()),
            "--output-dir",
            str(struct_dir),
            "--chr",
            args.chr,
            "--adaptive-qst",
            args.adaptive_qst,
            "--ve-ratios",
            args.ve_ratios,
            "--num-repeats",
            str(args.num_repeats),
            "--num-neutral",
            str(args.num_neutral),
            "--num-sim",
            str(args.num_sim),
            "--batch-size",
            str(args.batch_size),
            "--threshold-percentile",
            "0.95",
            "--num-pop",
            "2",
            "--num-ind",
            str(n_ind),
            "--num-rep",
            str(n_rep),
            "--code-dir",
            str(Path(args.code_dir).resolve()),
            "--fst-input-dir",
            str(Path(args.fst_input_dir).resolve()),
            "--r-env-tarball",
            args.r_env_tarball,
            "--sub-name",
            args.sub_name,
            "--request-memory",
            "4GB",
            "--mc-cores",
            "1",
        ]
        print(f"[{struct_id}] generating DAGs ...")
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
