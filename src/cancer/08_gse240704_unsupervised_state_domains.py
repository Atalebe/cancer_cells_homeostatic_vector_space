from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

from _gse240704_utils import load_config


def main(config_path: str) -> None:
    cfg = load_config(config_path)

    state_path = cfg["outputs"]["state_table"]
    out_dir = Path("results/gse240704/state_domains")
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(state_path)
    features = ["H", "S", "M", "R"]
    X = df[features].copy()

    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)

    rows = []
    best_k = None
    best_score = -np.inf
    best_labels = None

    for k in range(2, 9):
        model = KMeans(n_clusters=k, random_state=42, n_init=20)
        labels = model.fit_predict(Xs)
        score = silhouette_score(Xs, labels)
        rows.append({"k": k, "silhouette": float(score)})
        print(f"[info] k={k} silhouette={score:.6f}")
        if score > best_score:
            best_score = score
            best_k = k
            best_labels = labels

    sil_df = pd.DataFrame(rows)
    sil_df.to_csv(out_dir / "kmeans_silhouette_scan.csv", index=False)

    df["state_domain"] = [f"D{int(x)+1}" for x in best_labels]
    centroids = df.groupby("state_domain")[features].median().reset_index()
    counts = df["state_domain"].value_counts().rename_axis("state_domain").reset_index(name="n")

    df.to_parquet(out_dir / "state_domains.parquet", index=False)
    df.to_csv(out_dir / "state_domains.csv", index=False)
    centroids.to_csv(out_dir / "state_domain_centroids.csv", index=False)
    counts.to_csv(out_dir / "state_domain_counts.csv", index=False)

    print(f"[ok] wrote state domains to {out_dir}")
    print(f"[info] selected k: {best_k}")
    print(f"[info] best silhouette: {best_score:.6f}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    args = parser.parse_args()
    main(args.config)
