import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src.common.io import read_table, write_table
from src.common.score_utils import (
    expression_wide_to_cell_by_gene,
    attach_gene_symbols_to_cell_by_gene,
)


def safe_corr(x: pd.Series, y: pd.Series) -> float:
    x = pd.Series(x).astype(float)
    y = pd.Series(y).astype(float)
    if x.std(ddof=0) == 0 or y.std(ddof=0) == 0:
        return 0.0
    return float(x.corr(y))


def sector_vs_rest_logfc(x: pd.DataFrame, mask: pd.Series) -> pd.Series:
    """
    x is cell x gene log1p expression
    mask is boolean index aligned to rows of x
    """
    in_group = x.loc[mask]
    out_group = x.loc[~mask]

    mean_in = in_group.mean(axis=0)
    mean_out = out_group.mean(axis=0)

    # already in log1p space, so difference is a log-scale contrast
    return mean_in - mean_out


def build_top_table(metric: pd.Series, top_n: int = 50, ascending: bool = False) -> pd.DataFrame:
    s = metric.sort_values(ascending=ascending).head(top_n)
    return pd.DataFrame(
        {
            "gene": s.index,
            "score": s.values,
        }
    ).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--top-n", type=int, default=50)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    tables_dir = Path(cfg["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    expr_path = processed_dir / "expression_matrix_wide.parquet"
    gene_map_path = processed_dir / "gene_map.csv"
    refined_state_path = processed_dir / "state_table_refined.csv"
    branch_path = tables_dir / "branch_sector_assignments.csv"

    if not expr_path.exists():
        raise FileNotFoundError(f"Missing expression matrix: {expr_path}")
    if not gene_map_path.exists():
        raise FileNotFoundError(f"Missing gene map: {gene_map_path}")
    if not refined_state_path.exists():
        raise FileNotFoundError(f"Missing refined state table: {refined_state_path}")
    if not branch_path.exists():
        raise FileNotFoundError(f"Missing branch sector assignments: {branch_path}")

    expr = read_table(expr_path)
    gene_map = read_table(gene_map_path)
    state = read_table(refined_state_path)
    branch = read_table(branch_path)

    x = expression_wide_to_cell_by_gene(expr, gene_col="ensembl_id")
    x = attach_gene_symbols_to_cell_by_gene(
        x,
        gene_map,
        ensembl_col="ensembl_id",
        symbol_col="gene_symbol",
    )

    # align state/branch with expression matrix row order
    meta = state.merge(
        branch[["cell_id", "branch_sector"]],
        on="cell_id",
        how="left",
    ).set_index("cell_id")

    x = x.loc[meta.index]
    x_log = np.log1p(x)

    # ------------------------------------------------------------
    # 1. gene correlations with HRSM coordinates
    # ------------------------------------------------------------
    corr_rows = []
    for gene in x_log.columns:
        g = x_log[gene]
        corr_rows.append(
            {
                "gene": gene,
                "corr_H": safe_corr(g, meta["H"]),
                "corr_S": safe_corr(g, meta["S"]),
                "corr_M": safe_corr(g, meta["M"]),
                "corr_R": safe_corr(g, meta["R"]),
            }
        )

    corr_df = pd.DataFrame(corr_rows)
    corr_path = tables_dir / "gene_coordinate_correlations.csv"
    write_table(corr_df, corr_path)

    # top positive and negative per coordinate
    top_tables = {}
    for coord in ["H", "S", "M", "R"]:
        col = f"corr_{coord}"
        top_pos = corr_df.sort_values(col, ascending=False).head(args.top_n).copy()
        top_neg = corr_df.sort_values(col, ascending=True).head(args.top_n).copy()

        pos_path = tables_dir / f"top_positive_genes_corr_{coord}.csv"
        neg_path = tables_dir / f"top_negative_genes_corr_{coord}.csv"
        write_table(top_pos, pos_path)
        write_table(top_neg, neg_path)

        top_tables[f"top_positive_{coord}"] = str(pos_path)
        top_tables[f"top_negative_{coord}"] = str(neg_path)

    # ------------------------------------------------------------
    # 2. sector-vs-rest contrasts
    # ------------------------------------------------------------
    sector_tables = {}
    sectors_of_interest = [
        "malignant_reservoir",
        "unstable_committed_malignant",
        "stable_recoverable",
        "recoverable_reference",
        "transitional_intermediate",
    ]

    for sector in sectors_of_interest:
        mask = meta["branch_sector"] == sector
        n_in = int(mask.sum())
        n_out = int((~mask).sum())

        if n_in == 0 or n_out == 0:
            continue

        logfc = sector_vs_rest_logfc(x_log, mask)
        sector_df = pd.DataFrame(
            {
                "gene": logfc.index,
                "log1p_mean_diff_sector_vs_rest": logfc.values,
            }
        ).sort_values("log1p_mean_diff_sector_vs_rest", ascending=False)

        out_path = tables_dir / f"sector_vs_rest_{sector}.csv"
        write_table(sector_df, out_path)
        sector_tables[sector] = str(out_path)

        top_up_path = tables_dir / f"top_up_genes_{sector}.csv"
        top_down_path = tables_dir / f"top_down_genes_{sector}.csv"
        write_table(sector_df.head(args.top_n), top_up_path)
        write_table(sector_df.tail(args.top_n).sort_values("log1p_mean_diff_sector_vs_rest"), top_down_path)

        sector_tables[f"{sector}_top_up"] = str(top_up_path)
        sector_tables[f"{sector}_top_down"] = str(top_down_path)

    # ------------------------------------------------------------
    # 3. direct PKH26_High vs G1 contrast
    # ------------------------------------------------------------
    if {"PKH26_High", "G1"}.issubset(set(meta["population_label"].unique())):
        mask_high = meta["population_label"] == "PKH26_High"
        mask_g1 = meta["population_label"] == "G1"

        mean_high = x_log.loc[mask_high].mean(axis=0)
        mean_g1 = x_log.loc[mask_g1].mean(axis=0)
        diff = mean_high - mean_g1

        high_vs_g1 = pd.DataFrame(
            {
                "gene": diff.index,
                "log1p_mean_diff_PKH26_High_vs_G1": diff.values,
            }
        ).sort_values("log1p_mean_diff_PKH26_High_vs_G1", ascending=False)

        high_vs_g1_path = tables_dir / "PKH26_High_vs_G1_gene_contrast.csv"
        write_table(high_vs_g1, high_vs_g1_path)

        write_table(
            high_vs_g1.head(args.top_n),
            tables_dir / "top_up_genes_PKH26_High_vs_G1.csv",
        )
        write_table(
            high_vs_g1.tail(args.top_n).sort_values("log1p_mean_diff_PKH26_High_vs_G1"),
            tables_dir / "top_down_genes_PKH26_High_vs_G1.csv",
        )
    else:
        high_vs_g1_path = None

    # ------------------------------------------------------------
    # 4. summary
    # ------------------------------------------------------------
    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "n_cells": int(x_log.shape[0]),
        "n_genes": int(x_log.shape[1]),
        "outputs": {
            "gene_coordinate_correlations": str(corr_path),
            "top_coordinate_tables": top_tables,
            "sector_tables": sector_tables,
            "PKH26_High_vs_G1_gene_contrast": str(high_vs_g1_path) if high_vs_g1_path is not None else None,
        },
        "sector_counts": meta["branch_sector"].value_counts().to_dict(),
    }
    with open(tables_dir / "gene_program_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
