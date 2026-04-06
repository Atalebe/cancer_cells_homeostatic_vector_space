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
    score_gene_set_mean_z,
)


def safe_corr(x: pd.Series, y: pd.Series) -> float:
    x = pd.Series(x).astype(float)
    y = pd.Series(y).astype(float)
    if x.std(ddof=0) == 0 or y.std(ddof=0) == 0:
        return 0.0
    return float(x.corr(y))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--min-lists", type=int, default=3)
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    tables_dir = Path(cfg["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    expr = read_table(processed_dir / "expression_matrix_wide.parquet")
    gene_map = read_table(processed_dir / "gene_map.csv")
    state = read_table(processed_dir / "state_table_refined.csv")
    branch = read_table(tables_dir / "branch_sector_assignments.csv")
    core = read_table(tables_dir / "irreversible_core_candidate_top_ranked.csv")

    core_sub = core[core["n_lists_present"] >= args.min_lists].head(args.top_n).copy()
    selected_genes = core_sub["gene"].astype(str).tolist()

    x = expression_wide_to_cell_by_gene(expr, gene_col="ensembl_id")
    x = attach_gene_symbols_to_cell_by_gene(
        x,
        gene_map,
        ensembl_col="ensembl_id",
        symbol_col="gene_symbol",
    )

    core_score = score_gene_set_mean_z(x, selected_genes)

    merged = state.merge(branch[["cell_id", "branch_sector"]], on="cell_id", how="left")
    merged["candidate_core_score"] = merged["cell_id"].map(core_score.to_dict())

    out_path = tables_dir / "candidate_core_score_table.csv"
    write_table(merged, out_path)

    corr_summary = {
        "corr_core_H": safe_corr(merged["candidate_core_score"], merged["H"]),
        "corr_core_S": safe_corr(merged["candidate_core_score"], merged["S"]),
        "corr_core_M": safe_corr(merged["candidate_core_score"], merged["M"]),
        "corr_core_R": safe_corr(merged["candidate_core_score"], merged["R"]),
    }

    group_summary = (
        merged.groupby("population_label")[["candidate_core_score", "H", "S", "M", "R"]]
        .mean()
        .round(4)
        .reset_index()
    )
    group_summary_path = tables_dir / "candidate_core_score_group_summary.csv"
    write_table(group_summary, group_summary_path)

    sector_summary = (
        merged.groupby("branch_sector")[["candidate_core_score", "H", "S", "M", "R"]]
        .mean()
        .round(4)
        .reset_index()
        .sort_values("candidate_core_score", ascending=False)
    )
    sector_summary_path = tables_dir / "candidate_core_score_sector_summary.csv"
    write_table(sector_summary, sector_summary_path)

    corr_df = pd.DataFrame(
        [
            {
                "metric": k,
                "value": v,
            }
            for k, v in corr_summary.items()
        ]
    )
    corr_path = tables_dir / "candidate_core_score_correlations.csv"
    write_table(corr_df, corr_path)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "min_lists": int(args.min_lists),
        "top_n": int(args.top_n),
        "n_selected_genes": int(len(selected_genes)),
        "selected_genes": selected_genes,
        "correlations": corr_summary,
        "outputs": {
            "table": str(out_path),
            "group_summary": str(group_summary_path),
            "sector_summary": str(sector_summary_path),
            "correlations": str(corr_path),
        },
    }

    with open(tables_dir / "candidate_core_score_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
