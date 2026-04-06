import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

from src.common.io import read_table, write_table
from src.common.score_utils import (
    expression_wide_to_cell_by_gene,
    attach_gene_symbols_to_cell_by_gene,
    score_gene_set_mean_z,
)


def safe_corr(x: pd.Series, y: pd.Series) -> float:
    if x.std(ddof=0) == 0 or y.std(ddof=0) == 0:
        return 0.0
    return float(x.corr(y))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--min-lists", type=int, default=3)
    parser.add_argument("--top-n", type=int, default=16)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    tables_dir = Path(cfg["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    expr = read_table(processed_dir / "expression_matrix_wide.parquet")
    gene_map = read_table(processed_dir / "gene_map.csv")
    state = read_table(processed_dir / "state_table_refined.csv")
    core = read_table(tables_dir / "irreversible_core_candidate_top_ranked.csv")

    genes = core[core["n_lists_present"] >= args.min_lists].head(args.top_n)["gene"].astype(str).tolist()

    x = expression_wide_to_cell_by_gene(expr, gene_col="ensembl_id")
    x = attach_gene_symbols_to_cell_by_gene(
        x,
        gene_map,
        ensembl_col="ensembl_id",
        symbol_col="gene_symbol",
    )

    rows = []

    full_score = score_gene_set_mean_z(x, genes)
    full_df = state.copy()
    full_df["candidate_core_score"] = full_df["cell_id"].map(full_score.to_dict())

    rows.append(
        {
            "left_out_gene": "NONE",
            "n_genes_used": len(genes),
            "corr_H": safe_corr(full_df["candidate_core_score"], full_df["H"]),
            "corr_S": safe_corr(full_df["candidate_core_score"], full_df["S"]),
            "corr_M": safe_corr(full_df["candidate_core_score"], full_df["M"]),
            "corr_R": safe_corr(full_df["candidate_core_score"], full_df["R"]),
        }
    )

    for gene in genes:
        subgenes = [g for g in genes if g != gene]
        score = score_gene_set_mean_z(x, subgenes)
        df = state.copy()
        df["candidate_core_score"] = df["cell_id"].map(score.to_dict())

        rows.append(
            {
                "left_out_gene": gene,
                "n_genes_used": len(subgenes),
                "corr_H": safe_corr(df["candidate_core_score"], df["H"]),
                "corr_S": safe_corr(df["candidate_core_score"], df["S"]),
                "corr_M": safe_corr(df["candidate_core_score"], df["M"]),
                "corr_R": safe_corr(df["candidate_core_score"], df["R"]),
            }
        )

    out_df = pd.DataFrame(rows)
    out_path = tables_dir / "candidate_core_leave_one_gene_out.csv"
    write_table(out_df, out_path)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "top_n": args.top_n,
        "n_genes_tested": len(genes),
        "output_file": str(out_path),
    }

    with open(tables_dir / "candidate_core_leave_one_gene_out_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
