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


def sector_vs_rest_logfc(x: pd.DataFrame, mask: pd.Series) -> pd.Series:
    in_group = x.loc[mask]
    out_group = x.loc[~mask]
    return in_group.mean(axis=0) - out_group.mean(axis=0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--n-perm", type=int, default=200)
    parser.add_argument("--min-lists", type=int, default=3)
    parser.add_argument("--top-n", type=int, default=16)
    parser.add_argument("--sector", default="malignant_reservoir")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    tables_dir = Path(cfg["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    expr = read_table(processed_dir / "expression_matrix_wide.parquet")
    gene_map = read_table(processed_dir / "gene_map.csv")
    branch = read_table(tables_dir / "branch_sector_assignments.csv")
    core = read_table(tables_dir / "irreversible_core_candidate_top_ranked.csv")

    genes = core[core["n_lists_present"] >= args.min_lists].head(args.top_n)["gene"].astype(str).tolist()

    x = expression_wide_to_cell_by_gene(expr, gene_col="ensembl_id")
    x = attach_gene_symbols_to_cell_by_gene(
        x,
        gene_map,
        ensembl_col="ensembl_id",
        symbol_col="gene_symbol",
    )
    x_log = np.log1p(x)

    meta = branch.set_index("cell_id")
    x_log = x_log.loc[meta.index]

    real_mask = meta["branch_sector"] == args.sector
    real_logfc = sector_vs_rest_logfc(x_log[genes], real_mask)

    perm_rows = []
    labels = meta["branch_sector"].values.copy()

    for i in range(args.n_perm):
        shuffled = rng.permutation(labels)
        mask = pd.Series(shuffled == args.sector, index=meta.index)
        logfc = sector_vs_rest_logfc(x_log[genes], mask)
        for gene in genes:
            perm_rows.append(
                {
                    "perm_id": i,
                    "gene": gene,
                    "perm_logfc": float(logfc[gene]),
                }
            )

    perm_df = pd.DataFrame(perm_rows)

    summary_rows = []
    for gene in genes:
        null_vals = perm_df.loc[perm_df["gene"] == gene, "perm_logfc"]
        real_val = float(real_logfc[gene])
        z = 0.0
        if null_vals.std(ddof=0) > 0:
            z = (real_val - null_vals.mean()) / null_vals.std(ddof=0)
        summary_rows.append(
            {
                "gene": gene,
                "real_logfc": real_val,
                "null_mean": float(null_vals.mean()),
                "null_std": float(null_vals.std(ddof=0)),
                "z_score": float(z),
            }
        )

    summary_df = pd.DataFrame(summary_rows).sort_values("z_score", ascending=False)

    perm_path = tables_dir / f"sector_null_permutations_{args.sector}.csv"
    summary_path = tables_dir / f"sector_null_summary_{args.sector}.csv"
    write_table(perm_df, perm_path)
    write_table(summary_df, summary_path)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "sector": args.sector,
        "n_perm": args.n_perm,
        "top_n": args.top_n,
        "output_permutations": str(perm_path),
        "output_summary": str(summary_path),
    }

    with open(tables_dir / f"sector_null_summary_{args.sector}.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
