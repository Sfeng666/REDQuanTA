#!/usr/bin/env python3
"""Compare Stage C REDQuanTEA outputs to the reference production trees.

Usage:
  python scripts/compare_stageC.py

Exit 0 if Detection adaptive calls match and QST/TPR-FPR stay within tolerance.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

WORKSPACE = Path("/home/sfeng77/jobs/qst_workflow")
REF_DETECT = (
    WORKSPACE
    / "REDQuanTA_test_va_obs_softclip_floor_VGBVGW_flatdag_middlepeakremoved_correctva_abcVEVG"
    / "htcondor/results/compare_neutral_qst_strategy_floorwithoutsdenv"
    / "module2_original_bothnegna/investigate_anovenorep_hightpr"
    / "module1_F3_QST_ratioVbetweenVtotal/protein/summary/qst_results_abc.csv"
)
REF_DESIGN = (
    WORKSPACE
    / "REDQuanTA_test_va_obs_softclip_floor_VGBVGW_flatdag_middlepeakremoved_correctva_abcVEVG"
    / "htcondor/results/compare_neutral_qst_strategy_floorwithoutsdenv"
    / "module2_original_bothnegna/investigate_anovenorep_hightpr/production_F3"
)
NEW_DETECT = WORKSPACE / "REDQuanTEA/results/stageC_detect_protein/summary/qst_results_abc.csv"
NEW_DESIGN = WORKSPACE / "REDQuanTEA/results/stageC_design_production_F3"

QST_ATOL = 1e-4
TPR_ATOL = 0.01


def compare_detect() -> bool:
    if not NEW_DETECT.is_file():
        print(f"MISSING {NEW_DETECT}")
        return False
    ref = pd.read_csv(REF_DETECT)
    new = pd.read_csv(NEW_DETECT)
    ref["trait_id"] = ref["trait_id"].astype(str)
    new["trait_id"] = new["trait_id"].astype(str)
    m = ref.merge(new, on="trait_id", suffixes=("_ref", "_new"))
    print(f"Detection: {len(ref)} ref, {len(new)} new, {len(m)} merged")
    if len(m) != len(ref) or len(new) != len(ref):
        print("FAIL: trait set mismatch")
        only_ref = set(ref.trait_id) - set(new.trait_id)
        only_new = set(new.trait_id) - set(ref.trait_id)
        print("  only ref", len(only_ref), "only new", len(only_new))
        return False
    d_qst = (m["QST_new"] - m["QST_ref"]).abs()
    print(
        f"  |dQST| max={d_qst.max():.6g} median={d_qst.median():.6g} "
        f"n_gt_{QST_ATOL}={(d_qst > QST_ATOL).sum()}"
    )
    ad_ref = m["adaptive_ref"].astype(str).str.lower()
    ad_new = m["adaptive_new"].astype(str).str.lower()
    n_mismatch = int((ad_ref != ad_new).sum())
    print(f"  adaptive mismatches: {n_mismatch}")
    ok = n_mismatch == 0 and bool((d_qst <= QST_ATOL).all())
    print("Detection:", "PASS" if ok else "FAIL")
    return ok


def _load_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    ve = [c for c in df.columns if c.startswith("VEratio_")]
    return df.set_index("type")[ve]


def compare_design_chr(chr_name: str) -> bool:
    ref = REF_DESIGN / chr_name / f"tpr_fpr_matrix_{chr_name}_combo_0001.csv"
    new = NEW_DESIGN / chr_name / f"tpr_fpr_matrix_{chr_name}_combo_0001.csv"
    if not new.is_file():
        print(f"MISSING {new}")
        return False
    r = _load_matrix(ref)
    n = _load_matrix(new)
    common_rows = r.index.intersection(n.index)
    common_cols = r.columns.intersection(n.columns)
    if len(common_rows) == 0 or len(common_cols) == 0:
        print(f"FAIL {chr_name}: no overlapping matrix cells")
        return False
    d = (n.loc[common_rows, common_cols] - r.loc[common_rows, common_cols]).abs()
    print(
        f"Design {chr_name}: |d| max={d.max().max():.6g} "
        f"n_gt_{TPR_ATOL}={int((d > TPR_ATOL).sum().sum())} "
        f"shape ref={r.shape} new={n.shape}"
    )
    ok = bool((d <= TPR_ATOL).all().all())
    print(f"Design {chr_name}:", "PASS" if ok else "FAIL")
    return ok


def main() -> int:
    ok = True
    ok = compare_detect() and ok
    ok = compare_design_chr("autosomes") and ok
    ok = compare_design_chr("chrX") and ok
    print("STAGE C:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
