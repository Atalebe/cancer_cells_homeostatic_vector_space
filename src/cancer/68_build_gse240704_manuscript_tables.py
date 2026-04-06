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
    ap.add_argument("--outdir", default="results/gse240704/manuscript_tables")
    args = ap.parse_args()

    base = Path(args.base)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Table 1: geometry and domains
    geom = safe_read(base / "geometry" / "geometry_summary.csv")
    dom = safe_read(base / "state_domains" / "state_domain_counts.csv")
    if not geom.empty and not dom.empty:
        dom["table_section"] = "state_domains"
        dom.to_csv(outdir / "table1_state_domain_counts.csv", index=False)

    # Table 2: curated MGMT axis tests
    cond = safe_read(base / "condition_axis_tests_curated" / "condition_axis_tests_curated.csv")
    if not cond.empty:
        cond = cond.sort_values("p_value")
        cond.to_csv(outdir / "table2_curated_condition_axis_tests.csv", index=False)

    # Table 3: stable shell annotation
    stable = safe_read(base / "condition_axis_tests_curated" / "stable_shell_mgmt_overlay_repaired.csv")
    if not stable.empty:
        stable.to_csv(outdir / "table3_stable_shell_annotations.csv", index=False)

    # Table 4: domain enrichment
    enr = safe_read(base / "domain_condition_followup" / "domain_mgmt_enrichment_fisher.csv")
    if not enr.empty:
        enr = enr.sort_values("fisher_p_value")
        enr.to_csv(outdir / "table4_domain_mgmt_enrichment.csv", index=False)

    # Table 5: pairwise domain contrasts
    pw = safe_read(base / "domain_pairwise_followup" / "pairwise_domain_contrasts.csv")
    if not pw.empty:
        pw = pw.sort_values(["axis", "p_value"])
        pw.to_csv(outdir / "table5_pairwise_domain_contrasts.csv", index=False)

    # Table 6: D3 audit
    d3 = safe_read(base / "d3_state_audit" / "d3_metric_tests.csv")
    if not d3.empty:
        d3 = d3.sort_values("p_value")
        d3.to_csv(outdir / "table6_d3_metric_tests.csv", index=False)

    # Table 7: D3 shell proximity
    prox = safe_read(base / "d3_interpretation_followup" / "d3_shell_proximity_stratified_tests.csv")
    if not prox.empty:
        prox.to_csv(outdir / "table7_d3_shell_proximity_tests.csv", index=False)

    # Table 8: D3 directional chromosome arms
    arm = safe_read(base / "chromosome_arm_enrichment_followup" / "D3_directional_chromosome_arm_enrichment_combined.csv")
    if not arm.empty:
        arm = arm.sort_values(["direction", "fisher_p_value"])
        arm.to_csv(outdir / "table8_d3_chromosome_arm_enrichment.csv", index=False)

    # Table 9: top S-axis genes
    sgenes = safe_read(base / "final_directional_reports" / "s_axis_positive_loading_top_genes.csv")
    if not sgenes.empty:
        sgenes.to_csv(outdir / "table9_s_axis_positive_loading_top_genes.csv", index=False)

    # Table 10: D3 directional summary
    d3rep = safe_read(base / "final_directional_reports" / "D3_vs_not_D3_directional_compact_report.csv")
    if not d3rep.empty:
        d3rep.to_csv(outdir / "table10_d3_directional_compact_report.csv", index=False)

    print("[ok] wrote manuscript tables to", outdir)


if __name__ == "__main__":
    main()
