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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--k", type=int, default=10)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    tables_dir = Path(cfg["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    proxy_path = processed_dir / "hrsm_proxy_table.csv"
    if not proxy_path.exists():
        raise FileNotFoundError(f"Missing proxy table: {proxy_path}")

    df = read_table(proxy_path).copy()
    pc_cols = ["pc1", "pc2", "pc3", "pc4", "pc5"]
    X = df[pc_cols].values
    labels = df["population_label"].astype(str).values

    n_cells = X.shape[0]
    k = min(args.k, max(2, n_cells - 1))

    nbrs = NearestNeighbors(n_neighbors=k + 1, metric="euclidean")
    nbrs.fit(X)
    distances, indices = nbrs.kneighbors(X)

    same_label_purity = []
    local_entropy = []

    for i in range(n_cells):
        neigh_idx = indices[i, 1:]  # drop self
        neigh_labels = pd.Series(labels[neigh_idx])
        counts = neigh_labels.value_counts()
        purity = float(counts.get(labels[i], 0) / len(neigh_idx))
        entropy = shannon_entropy_from_counts(counts)
        same_label_purity.append(purity)
        local_entropy.append(entropy)

    df["knn_same_label_purity"] = same_label_purity
    df["knn_label_entropy"] = local_entropy

    # centroid margin: distance to own centroid minus nearest competing centroid
    centroids = (
        df.groupby("population_label")[pc_cols]
        .mean()
        .to_dict(orient="index")
    )

    own_dist = []
    nearest_other_dist = []
    centroid_margin = []

    for _, row in df.iterrows():
        own_label = row["population_label"]
        vec = row[pc_cols].values.astype(float)

        dists = {}
        for label, centroid_dict in centroids.items():
            centroid_vec = np.array([centroid_dict[c] for c in pc_cols], dtype=float)
            dists[label] = float(np.sqrt(((vec - centroid_vec) ** 2).sum()))

        d_own = dists[own_label]
        d_other = min(v for k2, v in dists.items() if k2 != own_label)

        own_dist.append(d_own)
        nearest_other_dist.append(d_other)
        centroid_margin.append(d_other - d_own)

    df["own_centroid_dist"] = own_dist
    df["nearest_other_centroid_dist"] = nearest_other_dist
    df["centroid_margin"] = centroid_margin

    # refined S: high purity, low entropy, high centroid margin
    df["S_refined_raw"] = (
        0.45 * robust_zscore(df["knn_same_label_purity"])
        + 0.25 * robust_zscore(-df["knn_label_entropy"])
        + 0.30 * robust_zscore(df["centroid_margin"])
    )
    df["S_refined"] = robust_zscore(df["S_refined_raw"])

    out_cols = [
        "cell_id",
        "population_label",
        "group_order",
        "pc1",
        "pc2",
        "pc3",
        "pc4",
        "pc5",
        "knn_same_label_purity",
        "knn_label_entropy",
        "own_centroid_dist",
        "nearest_other_centroid_dist",
        "centroid_margin",
        "S_refined_raw",
        "S_refined",
    ]
    out_df = df[out_cols].copy()

    out_path = tables_dir / "refined_structural_stability_table.csv"
    write_table(out_df, out_path)

    group_summary = (
        out_df.groupby("population_label")[
            ["knn_same_label_purity", "knn_label_entropy", "centroid_margin", "S_refined"]
        ]
        .agg(["mean", "median", "std"])
        .round(4)
    )
    group_summary.columns = ["_".join(col).strip("_") for col in group_summary.columns]
    group_summary = group_summary.reset_index()

    group_summary_path = tables_dir / "refined_structural_stability_group_summary.csv"
    write_table(group_summary, group_summary_path)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "k_neighbors": int(k),
        "outputs": {
            "table": str(out_path),
            "group_summary": str(group_summary_path),
        },
        "group_means": (
            out_df.groupby("population_label")[["S_refined"]]
            .mean()
            .round(4)
            .to_dict(orient="index")
        ),
    }

    with open(tables_dir / "refined_structural_stability_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
