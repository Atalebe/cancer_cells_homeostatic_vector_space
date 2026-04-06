#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def main() -> None:
    cfg = read_config()
    proc_dir = REPO_ROOT / cfg["processed_dir"]
    out_dir = REPO_ROOT / cfg["results_dir"] / "state_domains"
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(proc_dir / "state_table.parquet")
    X = df[["H", "S", "M", "R"]].values

    rows = []
    best_k = None
    best_score = -1
    best_labels = None

    for k in range(2, 9):
        km = KMeans(n_clusters=k, random_state=cfg["seed"], n_init=20)
        labels = km.fit_predict(X)
        score = silhouette_score(X, labels)
        rows.append({"k": k, "silhouette": float(score)})
        if score > best_score:
            best_score = score
            best_k = k
            best_labels = labels

    sil_df = pd.DataFrame(rows)
    sil_df.to_csv(out_dir / "silhouette_scan.csv", index=False)

    df["state_domain"] = [f"D{int(x)+1}" for x in best_labels]
    df.to_parquet(out_dir / "state_domains.parquet", index=False)
    df["state_domain"].value_counts().rename_axis("state_domain").reset_index(name="n").to_csv(
        out_dir / "state_domain_counts.csv", index=False
    )

    with open(out_dir / "state_domain_summary.json", "w", encoding="utf-8") as f:
        json.dump(
            {"best_k": int(best_k), "best_silhouette": float(best_score)},
            f,
            indent=2,
        )

    print("[ok] wrote state domain files")


if __name__ == "__main__":
    main()
