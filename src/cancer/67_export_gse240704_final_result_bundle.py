#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import pandas as pd


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def safe_read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="results/gse240704")
    ap.add_argument("--outdir", default="results/gse240704/final_result_bundle")
    args = ap.parse_args()

    base = Path(args.base)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summary = {}

    files = {
        "geometry_summary": base / "geometry" / "geometry_summary.csv",
        "state_domain_counts": base / "state_domains" / "state_domain_counts.csv",
        "condition_axis_tests_curated": base / "condition_axis_tests_curated" / "condition_axis_tests_curated.csv",
        "stable_shell_overlay": base / "condition_axis_tests_curated" / "stable_shell_mgmt_overlay_repaired.csv",
        "domain_mgmt_enrichment": base / "domain_condition_followup" / "domain_mgmt_enrichment_fisher.csv",
        "stable_shell_vs_rest": base / "domain_condition_followup" / "stable_shell_vs_rest_curated_tests.csv",
        "pairwise_domain_contrasts": base / "domain_pairwise_followup" / "pairwise_domain_contrasts.csv",
        "d3_metric_tests": base / "d3_state_audit" / "d3_metric_tests.csv",
        "d3_mgmt_enrichment": base / "d3_state_audit" / "d3_mgmt_enrichment.csv",
        "d3_shell_proximity": base / "d3_interpretation_followup" / "d3_shell_proximity_stratified_tests.csv",
        "s_axis_driver_genes": base / "biology_interpretation_followup" / "s_axis_driver_gene_token_counts.csv",
        "d3_gene_tokens": base / "biology_interpretation_followup" / "D3_vs_not_D3_gene_token_counts.csv",
        "d3_chr_selected_bg": base / "chromosome_enrichment_vs_selected_background" / "D3_directional_chromosome_enrichment_vs_selected_background.csv",
        "d3_arm_selected_bg": base / "chromosome_arm_enrichment_followup" / "D3_directional_chromosome_arm_enrichment_combined.csv",
        "d3_directional_report": base / "final_directional_reports" / "D3_vs_not_D3_directional_compact_report.csv",
        "s_axis_positive_report": base / "final_directional_reports" / "s_axis_positive_loading_top_genes.csv",
    }

    inventory_rows = []
    for label, path in files.items():
        exists = path.exists()
        n_rows = None
        if exists and path.suffix.lower() == ".csv":
            try:
                n_rows = len(pd.read_csv(path))
            except Exception:
                n_rows = None
        inventory_rows.append(
            {
                "label": label,
                "path": str(path),
                "exists": exists,
                "n_rows": n_rows,
            }
        )

    inv_df = pd.DataFrame(inventory_rows)
    inv_df.to_csv(outdir / "final_bundle_inventory.csv", index=False)

    geom = safe_read_csv(files["geometry_summary"])
    if not geom.empty:
        summary["geometry"] = geom.iloc[0].to_dict()

    states = safe_read_csv(files["state_domain_counts"])
    if not states.empty:
        summary["state_domain_counts"] = states.to_dict(orient="records")

    stable = safe_read_csv(files["stable_shell_overlay"])
    if not stable.empty:
        summary["stable_shell_samples"] = stable.to_dict(orient="records")

    cond = safe_read_csv(files["condition_axis_tests_curated"])
    if not cond.empty:
        summary["curated_condition_axis_tests"] = cond.to_dict(orient="records")

    chr_df = safe_read_csv(files["d3_chr_selected_bg"])
    if not chr_df.empty:
        top_chr = (
            chr_df.sort_values(["direction", "fisher_p_value", "fisher_odds_ratio"], ascending=[True, True, False])
            .groupby("direction", as_index=False)
            .head(5)
        )
        summary["top_directional_chromosome_hits"] = top_chr.to_dict(orient="records")

    arm_df = safe_read_csv(files["d3_arm_selected_bg"])
    if not arm_df.empty:
        top_arm = (
            arm_df.sort_values(["direction", "fisher_p_value", "fisher_odds_ratio"], ascending=[True, True, False])
            .groupby("direction", as_index=False)
            .head(5)
        )
        summary["top_directional_arm_hits"] = top_arm.to_dict(orient="records")

    with (outdir / "final_bundle_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote final bundle inventory:", outdir / "final_bundle_inventory.csv")
    print("[ok] wrote final bundle summary:", outdir / "final_bundle_summary.json")


if __name__ == "__main__":
    main()
