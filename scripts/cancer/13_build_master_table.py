import argparse
import json
from pathlib import Path

import pandas as pd
import yaml


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    raw_dir = Path(cfg["paths"]["raw_dir"])
    processed_dir = Path(cfg["paths"]["processed_dir"])
    metadata_dir = Path(cfg["paths"]["metadata_dir"])

    processed_dir.mkdir(parents=True, exist_ok=True)

    xlsx_path = raw_dir / "GSE124989_Reads.xlsx"
    meta_path = metadata_dir / "cell_metadata_registry.csv"

    if not xlsx_path.exists():
        raise FileNotFoundError(f"Missing workbook: {xlsx_path}")
    if not meta_path.exists():
        raise FileNotFoundError(f"Missing metadata registry: {meta_path}")

    expr = pd.read_excel(xlsx_path, sheet_name="Blad1")
    meta = pd.read_csv(meta_path)

    expr = expr.rename(columns={expr.columns[0]: "ensembl_id"})
    expr["ensembl_id"] = expr["ensembl_id"].astype(str)

    cell_columns = [c for c in expr.columns if c != "ensembl_id"]
    missing = sorted(set(meta["cell_id"]) - set(cell_columns))
    extra = sorted(set(cell_columns) - set(meta["cell_id"]))

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "n_genes": int(expr.shape[0]),
        "n_cell_columns": int(len(cell_columns)),
        "n_metadata_cells": int(meta.shape[0]),
        "missing_in_expression": len(missing),
        "extra_in_expression": len(extra),
    }

    expr_path = processed_dir / "expression_matrix_wide.parquet"
    meta_out_path = processed_dir / "cell_metadata_registry.parquet"

    expr.to_parquet(expr_path, index=False)
    meta.to_parquet(meta_out_path, index=False)

    with open(processed_dir / "master_table_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
