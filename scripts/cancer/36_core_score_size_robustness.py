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
    parser.add_argument("--sizes", nargs="+", type=int, default=[12, 16, 20])
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

    core = core[core["n_lists_present"] >= args.min_lists].copy()

    x = expression_wide_to_cell_by_gene(expr, gene_col="ensembl_id")
    x = attach_gene_symbols_to_cell_by_gene(
        x,
        gene_map,
        ensembl_col="ensembl_id",
        symbol_col="gene_symbol",
    )

    rows = []
    for n in args.sizes:
        genes = core.head(n)["gene"].astype(str).tolist()
        score = score_gene_set_mean_z(x, genes)
        score_map = score.to_dict()

        df = state.copy()
        df["candidate_core_score"] = df["cell_id"].map(score_map)

        rows.append(
            {
                "top_n": n,
                "n_genes_used": len(genes),
                "corr_H": safe_corr(df["candidate_core_score"], df["H"]),
                "corr_S": safe_corr(df["candidate_core_score"], df["S"]),
                "corr_M": safe_corr(df["candidate_core_score"], df["M"]),
                "corr_R": safe_corr(df["candidate_core_score"], df["R"]),
            }
        )

    out_df = pd.DataFrame(rows)
    out_path = tables_dir / "candidate_core_score_size_robustness.csv"
    write_table(out_df, out_path)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "sizes_tested": args.sizes,
        "output_file": str(out_path),
    }

    with open(tables_dir / "candidate_core_score_size_robustness_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
