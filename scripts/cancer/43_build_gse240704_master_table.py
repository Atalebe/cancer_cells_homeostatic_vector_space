import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

from src.common.io import read_table, write_table


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--matrix-kind", choices=["normalized", "non_normalized"], default="normalized")
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    metadata_dir = Path(cfg["paths"]["metadata_dir"])
    processed_dir.mkdir(parents=True, exist_ok=True)

    expr = read_table(processed_dir / f"gse240704_expression_{args.matrix_kind}.parquet")
    meta = read_table(metadata_dir / f"gse240704_metadata_registry_{args.matrix_kind}.csv")

    expr_cells = [c for c in expr.columns if c != "gene_id"]
    meta_cells = meta["cell_id"].astype(str).tolist()

    missing_in_expr = sorted(set(meta_cells) - set(expr_cells))
    extra_in_expr = sorted(set(expr_cells) - set(meta_cells))

    expr_t = expr.set_index("gene_id").T.reset_index().rename(columns={"index": "cell_id"})
    merged = meta.merge(expr_t, on="cell_id", how="inner")

    out_parquet = processed_dir / f"gse240704_master_cell_by_gene_{args.matrix_kind}.parquet"
    write_table(merged, out_parquet)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "matrix_kind": args.matrix_kind,
        "n_genes": int(expr.shape[0]),
        "n_expression_cells": int(len(expr_cells)),
        "n_metadata_cells": int(len(meta_cells)),
        "n_merged_cells": int(merged.shape[0]),
        "missing_in_expression": len(missing_in_expr),
        "extra_in_expression": len(extra_in_expr),
        "output_parquet": str(out_parquet),
    }

    with open(processed_dir / f"gse240704_master_table_{args.matrix_kind}_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
