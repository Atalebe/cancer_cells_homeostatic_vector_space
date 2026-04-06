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


def rank_of_gene(corr_df: pd.DataFrame, gene: str, col: str, ascending: bool) -> int | None:
    df = corr_df.sort_values(col, ascending=ascending).reset_index(drop=True)
    hits = df.index[df["gene"] == gene].tolist()
    return int(hits[0] + 1) if hits else None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--gene", default="ESYT2")
    args = parser.parse_args()

    gene = args.gene.upper().strip()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    tables_dir = Path(cfg["paths"]["tables_dir"])
    figures_dir = Path(cfg["paths"]["figures_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    expr = read_table(processed_dir / "expression_matrix_wide.parquet")
    gene_map = read_table(processed_dir / "gene_map.csv")
    state = read_table(processed_dir / "state_table_refined.csv")
    branch = read_table(tables_dir / "branch_sector_assignments.csv")
    corr_df = read_table(tables_dir / "gene_coordinate_correlations.csv")
    core_score = read_table(tables_dir / "candidate_core_score_table.csv")

    x = expression_wide_to_cell_by_gene(expr, gene_col="ensembl_id")
    x = attach_gene_symbols_to_cell_by_gene(
        x,
        gene_map,
        ensembl_col="ensembl_id",
        symbol_col="gene_symbol",
    )
    x_log = np.log1p(x)

    if gene not in {c.split("__dup")[0] for c in x_log.columns}:
        raise ValueError(f"{gene} not found in expression matrix after symbol mapping.")

    gene_cols = [c for c in x_log.columns if c.split("__dup")[0] == gene]
    if len(gene_cols) > 1:
        esyt = x_log[gene_cols].mean(axis=1)
    else:
        esyt = x_log[gene_cols[0]]

    meta = state.merge(branch[["cell_id", "branch_sector"]], on="cell_id", how="left")
    meta = meta.merge(core_score[["cell_id", "candidate_core_score"]], on="cell_id", how="left")
    meta["gene_log1p_expr"] = meta["cell_id"].map(esyt.to_dict())

    # summaries
    pop_summary = (
        meta.groupby("population_label")[["gene_log1p_expr", "H", "S", "M", "R", "candidate_core_score"]]
        .mean()
        .round(4)
        .reset_index()
    )
    pop_summary_path = tables_dir / f"{gene}_population_summary.csv"
    write_table(pop_summary, pop_summary_path)

    sector_summary = (
        meta.groupby("branch_sector")[["gene_log1p_expr", "H", "S", "M", "R", "candidate_core_score"]]
        .mean()
        .round(4)
        .reset_index()
        .sort_values("gene_log1p_expr", ascending=False)
    )
    sector_summary_path = tables_dir / f"{gene}_sector_summary.csv"
    write_table(sector_summary, sector_summary_path)

    correlations = {
        "corr_gene_H": safe_corr(meta["gene_log1p_expr"], meta["H"]),
        "corr_gene_S": safe_corr(meta["gene_log1p_expr"], meta["S"]),
        "corr_gene_M": safe_corr(meta["gene_log1p_expr"], meta["M"]),
        "corr_gene_R": safe_corr(meta["gene_log1p_expr"], meta["R"]),
        "corr_gene_core_score": safe_corr(meta["gene_log1p_expr"], meta["candidate_core_score"]),
    }
    corr_path = tables_dir / f"{gene}_correlations.csv"
    write_table(
        pd.DataFrame([{"metric": k, "value": v} for k, v in correlations.items()]),
        corr_path,
    )

    ranks = {
        "rank_positive_H": rank_of_gene(corr_df, gene, "corr_H", ascending=False),
        "rank_negative_S": rank_of_gene(corr_df, gene, "corr_S", ascending=True),
        "rank_positive_M": rank_of_gene(corr_df, gene, "corr_M", ascending=False),
        "rank_negative_R": rank_of_gene(corr_df, gene, "corr_R", ascending=True),
    }

    # top expressing cells
    top_cells = (
        meta[["cell_id", "population_label", "branch_sector", "gene_log1p_expr", "H", "S", "M", "R", "candidate_core_score"]]
        .sort_values("gene_log1p_expr", ascending=False)
        .head(25)
    )
    top_cells_path = tables_dir / f"{gene}_top_cells.csv"
    write_table(top_cells, top_cells_path)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "gene": gene,
        "gene_columns_used": gene_cols,
        "correlations": correlations,
        "ranks": ranks,
        "outputs": {
            "population_summary": str(pop_summary_path),
            "sector_summary": str(sector_summary_path),
            "correlations": str(corr_path),
            "top_cells": str(top_cells_path),
        },
    }

    with open(tables_dir / f"{gene}_diagnostic_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
