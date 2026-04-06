import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

from src.common.io import read_table, write_table
from src.common.qc_utils import (
    compute_cell_qc_from_wide_expression,
    summarize_qc_by_group,
    gene_id_duplicate_summary,
    summarize_missingness,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    metadata_dir = Path(cfg["paths"]["metadata_dir"])
    audit_dir = Path(cfg["paths"]["audit_dir"])
    audit_dir.mkdir(parents=True, exist_ok=True)

    expr_path = processed_dir / "expression_matrix_wide.parquet"
    meta_path = metadata_dir / "cell_metadata_registry.csv"

    if not expr_path.exists():
        raise FileNotFoundError(f"Missing processed expression matrix: {expr_path}")
    if not meta_path.exists():
        raise FileNotFoundError(f"Missing metadata registry: {meta_path}")

    expr = read_table(expr_path)
    meta = read_table(meta_path)

    qc = compute_cell_qc_from_wide_expression(expr, gene_col="ensembl_id")
    qc = qc.merge(meta, on="cell_id", how="left")

    qc_path = audit_dir / "cell_qc_metrics.csv"
    write_table(qc, qc_path)

    group_summary = summarize_qc_by_group(qc, group_col="population_label")
    group_summary_path = audit_dir / "group_qc_summary.csv"
    write_table(group_summary, group_summary_path)

    meta_missing = summarize_missingness(meta)
    meta_missing_path = audit_dir / "metadata_missingness.csv"
    write_table(meta_missing, meta_missing_path)

    dup_summary = gene_id_duplicate_summary(expr, gene_col="ensembl_id")

    audit_summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "n_cells": int(qc.shape[0]),
        "n_genes": int(expr.shape[0]),
        "population_counts": qc["population_label"].value_counts(dropna=False).to_dict(),
        "qc_outputs": {
            "cell_qc_metrics": str(qc_path),
            "group_qc_summary": str(group_summary_path),
            "metadata_missingness": str(meta_missing_path),
        },
        "total_counts_median": float(qc["total_counts"].median()),
        "detected_genes_median": float(qc["detected_genes"].median()),
        "zero_fraction_median": float(qc["zero_fraction"].median()),
        "gene_id_summary": dup_summary,
    }

    with open(audit_dir / "audit_summary.json", "w", encoding="utf-8") as fh:
        json.dump(audit_summary, fh, indent=2)

    print(json.dumps(audit_summary, indent=2))


if __name__ == "__main__":
    main()
