#!/usr/bin/env python3
"""Aggregate per-trait *_result.csv into qst_results_abc.csv and adaptive summary."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_OUT = SCRIPT_DIR.parent / "results" / "traits_harmonizr_validation"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--summary-dir",
        type=Path,
        default=None,
        help="Directory for qst_results_abc.csv and adaptive_trait_summary.csv (default: results-dir)",
    )
    args = parser.parse_args()
    summary_dir = args.summary_dir or args.results_dir
    summary_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for path in sorted(args.results_dir.glob("trait_*/*_result.csv")):
        tid = path.parent.name.replace("trait_", "", 1)
        with open(path, newline="") as f:
            r = next(csv.DictReader(f))
        r["trait_id"] = tid
        rows.append(r)

    if not rows:
        raise SystemExit(f"No *_result.csv under {args.results_dir}")

    out_csv = summary_dir / "qst_results_abc.csv"
    fields = list(rows[0].keys())
    if "trait_id" in fields:
        fields = ["trait_id"] + [f for f in fields if f != "trait_id"]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    n = len(rows)
    n_adaptive = sum(1 for r in rows if str(r.get("adaptive", "")).lower() == "yes")
    n_no = sum(1 for r in rows if str(r.get("adaptive", "")).lower() == "no")
    n_na = n - n_adaptive - n_no
    n_classifiable = n_adaptive + n_no

    summary = {
        "n_traits": n,
        "n_adaptive": n_adaptive,
        "n_non_adaptive": n_no,
        "n_na_adaptive": n_na,
        "n_classifiable": n_classifiable,
        "pct_adaptive": round(100 * n_adaptive / n_classifiable, 3) if n_classifiable else 0.0,
        "pct_non_adaptive": round(100 * n_no / n_classifiable, 3) if n_classifiable else 0.0,
    }
    summary_path = summary_dir / "adaptive_trait_summary.csv"
    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary.keys()))
        w.writeheader()
        w.writerow(summary)

    print(f"Wrote {out_csv} ({n} traits)")
    print(f"Adaptive: {n_adaptive} ({summary['pct_adaptive']}% of {n_classifiable} classifiable)")
    print(f"Non-adaptive: {n_no} ({summary['pct_non_adaptive']}% of classifiable)")
    if n_na:
        print(f"NA adaptive call: {n_na}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
