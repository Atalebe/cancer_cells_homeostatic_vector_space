#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import sparse


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def main() -> None:
    cfg = read_config()
    proc_dir = REPO_ROOT / cfg["processed_dir"]
    out_dir = REPO_ROOT / cfg["results_dir"] / "metadata_inspection"
    out_dir.mkdir(parents=True, exist_ok=True)

    mat = sparse.load_npz(proc_dir / "counts_sparse.npz")
    meta = pd.read_parquet(proc_dir / "cell_metadata_registry_aligned.parquet")

    total_counts = np.asarray(mat.sum(axis=0)).ravel()
    detected_genes = np.asarray((mat > 0).sum(axis=0)).ravel()

    summaries = {
        "n_genes": int(mat.shape[0]),
        "n_cells": int(mat.shape[1]),
        "n_nonzero_entries": int(mat.nnz),
        "median_total_counts": float(np.median(total_counts)),
        "median_detected_genes": float(np.median(detected_genes)),
        "n_metadata_rows": int(len(meta)),
        "n_cells_without_any_metadata_match": int((~meta.drop(columns=["cell_id", "gsm"]).notna().any(axis=1)).sum()),
        "metadata_columns": meta.columns.tolist(),
    }

    with open(out_dir / "realization_inspection_summary.json", "w", encoding="utf-8") as f:
        json.dump(summaries, f, indent=2)

    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
