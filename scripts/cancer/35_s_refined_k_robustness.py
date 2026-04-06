import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.neighbors import NearestNeighbors

from src.common.io import read_table, write_table
from src.common.qc_utils import robust_zscore


def shannon_entropy_from_counts(counts: pd.Series) -> float:
    probs = counts / counts.sum()
    probs = probs[probs > 0]
    return float(-(probs * np.log(probs)).sum())


def compute_s(df: pd.DataFrame, k: int) -> pd.DataFrame:
    pc_cols = ["pc1", "pc2", "pc3", "pc4", "pc5"]
    X = df[pc_cols].values
    labels = df["population_label"].astype(str).values

    n_cells = X.shape[0]
    k_eff = min(k, max(2, n_cells - 1))

    nbrs = NearestNeighbors(n_neighbors=k_eff + 1, metric="euclidean")
    nbrs.fit(X)
    _, indices = nbrs.kneighbors(X)

    purity_vals = []
    entropy_vals = []

    for i in range(n_cells):
        neigh_idx = indices[i, 1:]
        neigh_labels = pd.Series(labels[neigh_idx])
        counts = neigh_labels.value_counts()
        purity = float(counts.get(labels[i], 0) / len(neigh_idx))
        entropy = shannon_entropy_from_counts(counts)
        purity_vals.append(purity)
        entropy_vals.append(entropy)

    centroids = df.groupby("population_label")[pc_cols].mean().to_dict(orient="index")

    centroid_margin = []
    for _, row in df.iterrows():
        own_label = row["population_label"]
        vec = row[pc_cols].values.astype(float)

        dists = {}
        for label, centroid_dict in centroids.items():
            centroid_vec = np.array([centroid_dict[c] for c in pc_cols], dtype=float)
            dists[label] = float(np.sqrt(((vec - centroid_vec) ** 2).sum()))

        d_own = dists[own_label]
        d_other = min(v for kk, v in dists.items() if kk != own_label)
        centroid_margin.append(d_other - d_own)

    out = df[["cell_id", "population_label"]].copy()
    out["k"] = k_eff
    out["knn_same_label_purity"] = purity_vals
    out["knn_label_entropy"] = entropy_vals
    out["centroid_margin"] = centroid_margin
    out["S_refined_raw"] = (
        0.45 * robust_zscore(out["knn_same_label_purity"])
        + 0.25 * robust_zscore(-out["knn_label_entropy"])
        + 0.30 * robust_zscore(out["centroid_margin"])
    )
    out["S_refined"] = robust_zscore(out["S_refined_raw"])
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--k-values", nargs="+", type=int, default=[8, 10, 12, 15])
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    tables_dir = Path(cfg["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    proxy = read_table(processed_dir / "hrsm_proxy_table.csv")

    all_rows = []
    group_rows = []

    for k in args.k_values:
        out = compute_s(proxy, k)
        all_rows.append(out)

        grp = (
            out.groupby("population_label")[["S_refined", "knn_same_label_purity", "knn_label_entropy", "centroid_margin"]]
            .mean()
            .round(4)
            .reset_index()
        )
        grp["k"] = k
        group_rows.append(grp)

    full_df = pd.concat(all_rows, ignore_index=True)
    group_df = pd.concat(group_rows, ignore_index=True)

    full_path = tables_dir / "refined_s_k_robustness_full.csv"
    group_path = tables_dir / "refined_s_k_robustness_group_summary.csv"
    write_table(full_df, full_path)
    write_table(group_df, group_path)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "k_values": args.k_values,
        "output_full": str(full_path),
        "output_group_summary": str(group_path),
    }

    with open(tables_dir / "refined_s_k_robustness_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
