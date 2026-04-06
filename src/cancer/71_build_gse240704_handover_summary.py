#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import pandas as pd


def safe_read(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="results/gse240704")
    ap.add_argument("--outjson", default="results/gse240704/final_result_bundle/gse240704_handover_summary.json")
    args = ap.parse_args()

    base = Path(args.base)
    outjson = Path(args.outjson)
    outjson.parent.mkdir(parents=True, exist_ok=True)

    summary = {
        "dataset": "GSE240704",
        "branch_status": "analysis branch consolidated",
        "main_claims": [],
        "exact_pngs": [],
    }

    states = safe_read(base / "state_domains" / "state_domain_counts.csv")
    if not states.empty:
        summary["main_claims"].append({
            "claim": "Three unsupervised state domains were identified",
            "details": states.to_dict(orient="records"),
        })

    stable = safe_read(base / "condition_axis_tests_curated" / "stable_shell_mgmt_overlay_repaired.csv")
    if not stable.empty:
        summary["main_claims"].append({
            "claim": "Stable shell contains three samples, two MGMT-methylated and one unresolved",
            "details": stable[["sample_id", "state_domain", "placeholder_condition_cur", "phi", "H", "S", "M", "R"]].to_dict(orient="records"),
        })

    cond = safe_read(base / "condition_axis_tests_curated" / "condition_axis_tests_curated.csv")
    if not cond.empty:
        summary["main_claims"].append({
            "claim": "Curated MGMT comparison shows strongest shifts on H and R, weaker on phi, minimal on S and M",
            "details": cond.to_dict(orient="records"),
        })

    d3arm = safe_read(base / "chromosome_arm_enrichment_followup" / "D3_directional_chromosome_arm_enrichment_combined.csv")
    if not d3arm.empty:
        top = d3arm.sort_values(["direction", "fisher_p_value"]).groupby("direction", as_index=False).head(5)
        summary["main_claims"].append({
            "claim": "D3 has directional genomic structure, especially lower-in-D3 enrichment on 21q and higher-in-D3 enrichment on 7p and 4q",
            "details": top.to_dict(orient="records"),
        })

    sgenes = safe_read(base / "final_directional_reports" / "s_axis_positive_loading_top_genes.csv")
    if not sgenes.empty:
        summary["main_claims"].append({
            "claim": "Positive S-axis loading is driven by a compact annotated gene set led by PTCH1, NR3C1, PRDM16, and TP73",
            "details": sgenes.head(15).to_dict(orient="records"),
        })

    summary["exact_pngs"] = [
        "hm_plane.png",
        "hs_plane.png",
        "mr_plane.png",
        "phi_r_plane.png",
        "sr_plane.png",
        "mgmt_counts_by_state_domain.png",
        "h_by_mgmt_condition.png",
        "s_by_mgmt_condition.png",
        "m_by_mgmt_condition.png",
        "r_by_mgmt_condition.png",
        "phi_by_mgmt_condition.png",
        "h_by_state_domain.png",
        "s_by_state_domain.png",
        "m_by_state_domain.png",
        "r_by_state_domain.png",
        "phi_by_state_domain.png",
        "S_vs_std_beta_by_domain.png",
        "mean_beta_vs_std_beta_by_domain.png",
        "phi_vs_std_beta_by_domain.png",
        "iqr_beta_by_state_domain.png",
        "min_hrsm_distance_to_shell_by_state_domain.png",
        "std_beta_by_state_domain.png",
        "mean_beta_by_state_domain.png",
        "phi_vs_min_hrsm_distance_to_shell.png",
        "D3_directional_chromosome_log2OR_vs_selected_background.png",
        "D3_directional_chromosome_neglog10p_vs_selected_background.png",
        "D3_chromosome21_directional_log2OR.png",
    ]

    with outjson.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote handover summary:", outjson)


if __name__ == "__main__":
    main()
