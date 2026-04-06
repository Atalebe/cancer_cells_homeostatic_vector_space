import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from src.common.io import read_table, write_table
from src.common.qc_utils import robust_zscore


def inverse_centroid_distance_score(df: pd.DataFrame, label_col: str, pc_cols: list[str]) -> pd.Series:
    out = pd.Series(index=df.index, dtype=float)
    for label, sub in df.groupby(label_col):
        centroid = sub[pc_cols].mean(axis=0).values
        d = np.sqrt(((sub[pc_cols].values - centroid) ** 2).sum(axis=1))
        out.loc[sub.index] = -d
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    audit_dir = Path(cfg["paths"]["audit_dir"])
    results_tables = Path(cfg["paths"]["tables_dir"])
    results_tables.mkdir(parents=True, exist_ok=True)

    expr_path = processed_dir / "expression_matrix_wide.parquet"
    meta_path = processed_dir / "cell_metadata_registry.parquet"
    qc_path = audit_dir / "cell_qc_metrics.csv"

    expr = read_table(expr_path)
    meta = read_table(meta_path)
    qc = read_table(qc_path)

    cell_cols = [c for c in expr.columns if c != "ensembl_id"]
    X = expr[cell_cols].T
    X.index = cell_cols

    # mild filtering for PCA stability
    gene_nonzero_frac = (X > 0).mean(axis=0)
    keep_genes = gene_nonzero_frac[gene_nonzero_frac >= 0.05].index
    Xf = X[keep_genes]

    Xlog = np.log1p(Xf)
    Xscaled = StandardScaler(with_mean=True, with_std=True).fit_transform(Xlog)

    pca = PCA(n_components=5, random_state=0)
    pcs = pca.fit_transform(Xscaled)

    pca_df = pd.DataFrame(
        pcs,
        index=X.index,
        columns=[f"pc{i+1}" for i in range(pcs.shape[1])],
    ).reset_index().rename(columns={"index": "cell_id"})

    df = meta.merge(qc, on=["cell_id", "raw_prefix", "population_label", "group_order", "dataset_id", "accession"], how="left")
    df = df.merge(pca_df, on="cell_id", how="left")

    # H reserve: group order plus library complexity contribution
    df["H_raw"] = (
        0.65 * robust_zscore(df["group_order"])
        + 0.20 * robust_zscore(df["log_total_counts"])
        + 0.15 * robust_zscore(df["log_detected_genes"])
    )

    # M commitment: conservative biological ordering
    df["M_raw"] = robust_zscore(df["group_order"])

    # S structural coherence: inverse distance to group centroid in first two PCs
    df["S_raw"] = inverse_centroid_distance_score(df, "population_label", ["pc1", "pc2"])
    df["S_raw"] = robust_zscore(df["S_raw"])

    # R recoverability: anchored to lower-order baseline and penalized by deviation
    g1_pc1 = float(df.loc[df["population_label"] == "G1", "pc1"].median())
    g1_pc2 = float(df.loc[df["population_label"] == "G1", "pc2"].median())
    baseline_dist = np.sqrt((df["pc1"] - g1_pc1) ** 2 + (df["pc2"] - g1_pc2) ** 2)
    df["R_raw"] = robust_zscore(-baseline_dist) - 0.25 * robust_zscore(df["group_order"])

    proxy_cols = [
        "cell_id",
        "population_label",
        "group_order",
        "total_counts",
        "detected_genes",
        "zero_fraction",
        "log_total_counts",
        "log_detected_genes",
        "pc1",
        "pc2",
        "pc3",
        "pc4",
        "pc5",
        "H_raw",
        "S_raw",
        "M_raw",
        "R_raw",
    ]
    proxy_df = df[proxy_cols].copy()

    proxy_path = processed_dir / "hrsm_proxy_table.csv"
    write_table(proxy_df, proxy_path)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "n_cells": int(proxy_df.shape[0]),
        "n_pca_genes_used": int(len(keep_genes)),
        "explained_variance_ratio": {
            f"pc{i+1}": float(v) for i, v in enumerate(pca.explained_variance_ratio_)
        },
        "output_file": str(proxy_path),
    }

    with open(processed_dir / "hrsm_proxy_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
