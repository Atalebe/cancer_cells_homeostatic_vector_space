import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

from src.common.io import read_table, write_table


def merge_score_column(
    base: pd.DataFrame,
    path: Path,
    score_col: str,
    new_col: str,
) -> pd.DataFrame:
    df = read_table(path).copy()
    if "gene" not in df.columns or score_col not in df.columns:
        raise ValueError(f"Missing expected columns in {path}")
    df = df[["gene", score_col]].rename(
        columns={
            "gene": "gene_symbol",
            score_col: new_col,
        }
    )
    return base.merge(df, on="gene_symbol", how="left")


def interpret(row):
    tags = []
    if pd.notna(row.get("corr_H")) and row["corr_H"] > 0.4:
        tags.append("high_H")
    if pd.notna(row.get("corr_M")) and row["corr_M"] > 0.3:
        tags.append("high_M")
    if pd.notna(row.get("corr_R")) and row["corr_R"] < -0.3:
        tags.append("low_R")
    if pd.notna(row.get("malignant_reservoir_logfc")) and row["malignant_reservoir_logfc"] > 0.5:
        tags.append("reservoir_enriched")
    if pd.notna(row.get("unstable_committed_logfc")) and row["unstable_committed_logfc"] > 0.5:
        tags.append("unstable_committed_enriched")
    if pd.notna(row.get("stable_recoverable_logfc")) and row["stable_recoverable_logfc"] < -0.5:
        tags.append("stable_recoverable_depleted")
    return ";".join(tags)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--min-lists", type=int, default=2)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    tables_dir = Path(cfg["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    overlap_path = tables_dir / "irreversible_core_candidate_top_ranked.csv"
    corr_path = tables_dir / "gene_coordinate_correlations.csv"
    reservoir_path = tables_dir / "sector_vs_rest_malignant_reservoir.csv"
    unstable_path = tables_dir / "sector_vs_rest_unstable_committed_malignant.csv"
    stable_path = tables_dir / "sector_vs_rest_stable_recoverable.csv"
    high_vs_g1_path = tables_dir / "PKH26_High_vs_G1_gene_contrast.csv"

    if not overlap_path.exists():
        raise FileNotFoundError(f"Missing overlap table: {overlap_path}")

    base = read_table(overlap_path).copy()
    base = base[base["n_lists_present"] >= args.min_lists].copy()
    base = base.rename(columns={"gene": "gene_symbol"})

    if corr_path.exists():
        corr = read_table(corr_path).copy().rename(columns={"gene": "gene_symbol"})
        base = base.merge(corr, on="gene_symbol", how="left")

    if reservoir_path.exists():
        base = merge_score_column(
            base,
            reservoir_path,
            score_col="log1p_mean_diff_sector_vs_rest",
            new_col="malignant_reservoir_logfc",
        )

    if unstable_path.exists():
        base = merge_score_column(
            base,
            unstable_path,
            score_col="log1p_mean_diff_sector_vs_rest",
            new_col="unstable_committed_logfc",
        )

    if stable_path.exists():
        base = merge_score_column(
            base,
            stable_path,
            score_col="log1p_mean_diff_sector_vs_rest",
            new_col="stable_recoverable_logfc",
        )

    if high_vs_g1_path.exists():
        df = read_table(high_vs_g1_path).copy()
        if "gene" not in df.columns or "log1p_mean_diff_PKH26_High_vs_G1" not in df.columns:
            raise ValueError(f"Missing expected columns in {high_vs_g1_path}")
        df = df[["gene", "log1p_mean_diff_PKH26_High_vs_G1"]].rename(
            columns={
                "gene": "gene_symbol",
                "log1p_mean_diff_PKH26_High_vs_G1": "PKH26_High_vs_G1_logfc",
            }
        )
        base = base.merge(df, on="gene_symbol", how="left")

    base["interpretation_tags"] = base.apply(interpret, axis=1)

    preferred_cols = [
        "gene_symbol",
        "n_lists_present",
        "present_in_lists",
        "mean_rank_across_present_lists",
        "priority_score",
        "corr_H",
        "corr_S",
        "corr_M",
        "corr_R",
        "malignant_reservoir_logfc",
        "unstable_committed_logfc",
        "stable_recoverable_logfc",
        "PKH26_High_vs_G1_logfc",
        "interpretation_tags",
    ]
    out_df = base[preferred_cols].copy().sort_values(
        ["n_lists_present", "priority_score"],
        ascending=[False, False],
    )

    out_path = tables_dir / "candidate_core_interpretation_table.csv"
    write_table(out_df, out_path)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "min_lists": int(args.min_lists),
        "n_rows": int(out_df.shape[0]),
        "output_file": str(out_path),
    }

    with open(tables_dir / "candidate_core_interpretation_table_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
