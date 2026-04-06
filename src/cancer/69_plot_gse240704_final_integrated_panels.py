#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def safe_read(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def bar_plot(df: pd.DataFrame, x: str, y: str, outfile: Path, title: str, rotate: int = 45) -> None:
    plt.figure(figsize=(10, 5))
    plt.bar(df[x].astype(str), df[y])
    plt.title(title)
    plt.xticks(rotation=rotate, ha="right")
    plt.tight_layout()
    plt.savefig(outfile, dpi=180)
    plt.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="results/gse240704")
    ap.add_argument("--outdir", default="results/gse240704/final_summary_plots")
    args = ap.parse_args()

    base = Path(args.base)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    dom = safe_read(base / "state_domains" / "state_domain_counts.csv")
    if not dom.empty:
        bar_plot(dom, "state_domain", "n", outdir / "final_state_domain_counts.png", "State domain counts")

    cond = safe_read(base / "condition_axis_tests_curated" / "condition_axis_tests_curated.csv")
    if not cond.empty:
        use = cond[["axis", "delta_median_a_minus_b"]].copy()
        bar_plot(use, "axis", "delta_median_a_minus_b", outdir / "final_mgmt_axis_delta_medians.png", "MGMT methylated minus non-methylated median shifts", rotate=0)

    shell = safe_read(base / "domain_condition_followup" / "stable_shell_vs_rest_curated_tests.csv")
    if not shell.empty:
        axdf = shell[shell["comparison"].str.contains("stable_shell_vs_rest_axis_", na=False)].copy()
        if not axdf.empty:
            axdf["axis"] = axdf["comparison"].str.replace("stable_shell_vs_rest_axis_", "", regex=False)
            bar_plot(axdf, "axis", "delta_median_shell_minus_rest", outdir / "final_stable_shell_vs_rest_axis_deltas.png", "Stable shell versus rest axis shifts", rotate=0)

    chrdf = safe_read(base / "chromosome_enrichment_vs_selected_background" / "D3_directional_chromosome_enrichment_vs_selected_background.csv")
    if not chrdf.empty:
        higher = chrdf[chrdf["direction"] == "higher_in_a"].sort_values("fisher_p_value").head(10)
        lower = chrdf[chrdf["direction"] == "lower_in_a"].sort_values("fisher_p_value").head(10)

        if not higher.empty:
            higher = higher.assign(label=higher["chromosome"].astype(str))
            bar_plot(higher, "label", "fisher_odds_ratio", outdir / "final_d3_higher_chr_odds_ratio_top10.png", "Top chromosome odds ratios, D3 higher-in-a", rotate=0)
        if not lower.empty:
            lower = lower.assign(label=lower["chromosome"].astype(str))
            bar_plot(lower, "label", "fisher_odds_ratio", outdir / "final_d3_lower_chr_odds_ratio_top10.png", "Top chromosome odds ratios, D3 lower-in-a", rotate=0)

    arm = safe_read(base / "chromosome_arm_enrichment_followup" / "D3_directional_chromosome_arm_enrichment_combined.csv")
    if not arm.empty:
        top = arm.sort_values(["direction", "fisher_p_value"]).groupby("direction", as_index=False).head(8).copy()
        top["label"] = top["direction"] + ":" + top["chromosome_arm"]
        bar_plot(top, "label", "fisher_odds_ratio", outdir / "final_d3_arm_enrichment_top_hits.png", "Top chromosome-arm odds ratios by direction")

    print("[ok] wrote final summary plots to", outdir)


if __name__ == "__main__":
    main()
