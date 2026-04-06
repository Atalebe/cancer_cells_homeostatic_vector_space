#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def safe_read(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="results/gse240704")
    ap.add_argument("--outcsv", default="results/gse240704/final_result_bundle/gse240704_candidate_headline_table.csv")
    args = ap.parse_args()

    base = Path(args.base)
    outcsv = Path(args.outcsv)
    outcsv.parent.mkdir(parents=True, exist_ok=True)

    rows = []

    cond = safe_read(base / "condition_axis_tests_curated" / "condition_axis_tests_curated.csv")
    if not cond.empty:
        for _, r in cond.iterrows():
            rows.append({
                "headline_family": "MGMT_axis_shift",
                "label": str(r["axis"]),
                "effect": r.get("delta_median_a_minus_b"),
                "p_value": r.get("p_value"),
                "note": f"{r.get('condition_a')} minus {r.get('condition_b')}",
            })

    shell = safe_read(base / "domain_condition_followup" / "stable_shell_vs_rest_curated_tests.csv")
    if not shell.empty:
        use = shell[shell["comparison"].astype(str).str.contains("stable_shell_vs_rest_axis_", na=False)]
        for _, r in use.iterrows():
            rows.append({
                "headline_family": "Stable_shell_vs_rest",
                "label": str(r["comparison"]).replace("stable_shell_vs_rest_axis_", ""),
                "effect": r.get("delta_median_shell_minus_rest"),
                "p_value": r.get("p_value"),
                "note": "curated subset",
            })

    d3 = safe_read(base / "d3_interpretation_followup" / "d3_shell_proximity_stratified_tests.csv")
    if not d3.empty:
        use = d3[d3["comparison"].astype(str).str.contains("D3_shell_near_vs_far_", na=False)]
        for _, r in use.iterrows():
            rows.append({
                "headline_family": "D3_near_vs_far",
                "label": str(r["axis"]),
                "effect": r.get("delta_median_near_minus_far"),
                "p_value": r.get("p_value"),
                "note": "within D3 only",
            })

    arm = safe_read(base / "chromosome_arm_enrichment_followup" / "D3_directional_chromosome_arm_enrichment_combined.csv")
    if not arm.empty:
        top = arm.sort_values(["direction", "fisher_p_value"]).groupby("direction", as_index=False).head(5)
        for _, r in top.iterrows():
            rows.append({
                "headline_family": "D3_directional_arm_enrichment",
                "label": f"{r.get('direction')}:{r.get('chromosome_arm')}",
                "effect": r.get("fisher_odds_ratio"),
                "p_value": r.get("fisher_p_value"),
                "note": "selected-probe background",
            })

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["headline_family", "p_value"], na_position="last")
    out.to_csv(outcsv, index=False)
    print("[ok] wrote candidate headline table:", outcsv)


if __name__ == "__main__":
    main()
