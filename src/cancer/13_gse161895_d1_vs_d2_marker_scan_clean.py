#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def is_ercc(gene_id: str) -> bool:
    return str(gene_id).startswith("ERCC-")


def main() -> None:
    cfg = read_config()
    in_dir = REPO_ROOT / cfg["results_dir"] / "d1_d2_marker_scan"
    out_dir = REPO_ROOT / cfg["results_dir"] / "d1_d2_marker_scan_clean"
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(in_dir / "d1_vs_d2_marker_scan_full.csv")
    df["is_ercc"] = df["gene_id"].astype(str).map(is_ercc)

    clean = df.loc[~df["is_ercc"]].copy()
    clean = clean.sort_values(["q_value_bh", "abs_delta_mean"], ascending=[True, False]).reset_index(drop=True)

    d1_high = clean.sort_values(["delta_mean_d1_minus_d2", "q_value_bh"], ascending=[False, True]).head(200)
    d1_low = clean.sort_values(["delta_mean_d1_minus_d2", "q_value_bh"], ascending=[True, True]).head(200)

    clean.to_csv(out_dir / "d1_vs_d2_marker_scan_full_clean.csv", index=False)
    d1_high.to_csv(out_dir / "d1_high_top200_clean.csv", index=False)
    d1_low.to_csv(out_dir / "d1_low_top200_clean.csv", index=False)

    summary = {
        "n_genes_input": int(len(df)),
        "n_ercc_removed": int(df["is_ercc"].sum()),
        "n_genes_clean": int(len(clean)),
        "top_d1_high_first10": d1_high["gene_id"].head(10).tolist(),
        "top_d1_low_first10": d1_low["gene_id"].head(10).tolist(),
    }

    with open(out_dir / "d1_vs_d2_marker_scan_clean_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))
    print("\n[D1-high clean top 20]")
    print(d1_high.head(20)[["gene_id", "delta_mean_d1_minus_d2", "q_value_bh"]])
    print("\n[D1-low clean top 20]")
    print(d1_low.head(20)[["gene_id", "delta_mean_d1_minus_d2", "q_value_bh"]])


if __name__ == "__main__":
    main()
