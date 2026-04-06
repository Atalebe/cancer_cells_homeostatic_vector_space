from __future__ import annotations

import numpy as np
import pandas as pd


def summarize_missingness(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "column": df.columns,
            "n_missing": df.isna().sum().values,
            "fraction_missing": df.isna().mean().values,
        }
    ).sort_values(["fraction_missing", "n_missing"], ascending=False)


def robust_zscore(values: pd.Series) -> pd.Series:
    x = pd.Series(values, copy=True).astype(float)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    if mad == 0 or np.isnan(mad):
        return pd.Series(np.zeros(len(x)), index=x.index)
    return (x - med) / (1.4826 * mad)


def compute_cell_qc_from_wide_expression(
    expr: pd.DataFrame,
    gene_col: str = "ensembl_id",
) -> pd.DataFrame:
    cell_cols = [c for c in expr.columns if c != gene_col]
    mat = expr[cell_cols]

    total_counts = mat.sum(axis=0)
    detected_genes = (mat > 0).sum(axis=0)
    zero_fraction = (mat == 0).mean(axis=0)

    out = pd.DataFrame(
        {
            "cell_id": cell_cols,
            "total_counts": total_counts.values,
            "detected_genes": detected_genes.values,
            "zero_fraction": zero_fraction.values,
        }
    )
    out["log_total_counts"] = np.log1p(out["total_counts"])
    out["log_detected_genes"] = np.log1p(out["detected_genes"])
    return out


def summarize_qc_by_group(
    qc_df: pd.DataFrame,
    group_col: str = "population_label",
) -> pd.DataFrame:
    rows = []
    for group, sub in qc_df.groupby(group_col, dropna=False):
        rows.append(
            {
                group_col: group,
                "n_cells": int(len(sub)),
                "total_counts_median": float(sub["total_counts"].median()),
                "total_counts_mean": float(sub["total_counts"].mean()),
                "detected_genes_median": float(sub["detected_genes"].median()),
                "detected_genes_mean": float(sub["detected_genes"].mean()),
                "zero_fraction_median": float(sub["zero_fraction"].median()),
                "zero_fraction_mean": float(sub["zero_fraction"].mean()),
            }
        )
    return pd.DataFrame(rows).sort_values(group_col).reset_index(drop=True)


def gene_id_duplicate_summary(expr: pd.DataFrame, gene_col: str = "ensembl_id") -> dict:
    n_total = int(expr.shape[0])
    n_unique = int(expr[gene_col].astype(str).nunique())
    n_duplicated_rows = int(n_total - n_unique)
    return {
        "n_gene_rows": n_total,
        "n_unique_gene_ids": n_unique,
        "n_duplicated_gene_rows": n_duplicated_rows,
    }
