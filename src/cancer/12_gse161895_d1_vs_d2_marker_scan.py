#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from scipy.stats import mannwhitneyu


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = np.empty(n, dtype=float)

    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        q[i] = prev

    out = np.empty(n, dtype=float)
    out[order] = np.minimum(q, 1.0)
    return out


def main() -> None:
    cfg = read_config()
    proc_dir = REPO_ROOT / cfg["processed_dir"]
    state_dir = REPO_ROOT / cfg["results_dir"] / "state_domains"
    out_dir = REPO_ROOT / cfg["results_dir"] / "d1_d2_marker_scan"
    out_dir.mkdir(parents=True, exist_ok=True)

    mat = sparse.load_npz(proc_dir / "counts_sparse.npz").tocsc()
    gene_ids = pd.read_csv(proc_dir / "gene_ids.csv")["gene_id"].astype(str).tolist()
    cell_ids = pd.read_csv(proc_dir / "cell_ids.csv")["cell_id"].astype(str).tolist()

    state = pd.read_parquet(state_dir / "state_domains.parquet")[["cell_id", "state_domain"]]
    state = state.set_index("cell_id").reindex(cell_ids).reset_index()

    d1_idx = np.where(state["state_domain"].values == "D1")[0]
    d2_idx = np.where(state["state_domain"].values == "D2")[0]

    mat_d1 = mat[:, d1_idx]
    mat_d2 = mat[:, d2_idx]

    mean_d1 = np.asarray(mat_d1.mean(axis=1)).ravel()
    mean_d2 = np.asarray(mat_d2.mean(axis=1)).ravel()
    det_d1 = np.asarray((mat_d1 > 0).mean(axis=1)).ravel()
    det_d2 = np.asarray((mat_d2 > 0).mean(axis=1)).ravel()

    # restrict formal testing to moderately detected genes to save time
    keep = (det_d1 >= 0.05) | (det_d2 >= 0.05)
    keep_idx = np.where(keep)[0]

    rows = []
    pvals = []

    for i in keep_idx:
        a = mat_d1[i, :].toarray().ravel()
        b = mat_d2[i, :].toarray().ravel()

        try:
            u, p = mannwhitneyu(a, b, alternative="two-sided")
        except Exception:
            u, p = np.nan, 1.0

        rows.append({
            "gene_id": gene_ids[i],
            "mean_d1": float(mean_d1[i]),
            "mean_d2": float(mean_d2[i]),
            "delta_mean_d1_minus_d2": float(mean_d1[i] - mean_d2[i]),
            "detected_frac_d1": float(det_d1[i]),
            "detected_frac_d2": float(det_d2[i]),
            "delta_detected_frac_d1_minus_d2": float(det_d1[i] - det_d2[i]),
            "u_statistic": float(u) if not np.isnan(u) else None,
            "p_value": float(p),
        })
        pvals.append(p)

    out = pd.DataFrame(rows)
    out["q_value_bh"] = bh_adjust(np.array(pvals))
    out["abs_delta_mean"] = out["delta_mean_d1_minus_d2"].abs()

    out = out.sort_values(["q_value_bh", "abs_delta_mean"], ascending=[True, False]).reset_index(drop=True)
    out.to_csv(out_dir / "d1_vs_d2_marker_scan_full.csv", index=False)

    top_up = out.sort_values(["delta_mean_d1_minus_d2", "q_value_bh"], ascending=[False, True]).head(200)
    top_down = out.sort_values(["delta_mean_d1_minus_d2", "q_value_bh"], ascending=[True, True]).head(200)

    top_up.to_csv(out_dir / "d1_high_top200.csv", index=False)
    top_down.to_csv(out_dir / "d1_low_top200.csv", index=False)

    summary = {
        "n_genes_tested": int(len(out)),
        "n_d1_cells": int(len(d1_idx)),
        "n_d2_cells": int(len(d2_idx)),
    }
    with open(out_dir / "d1_vs_d2_marker_scan_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))
    print(top_up.head(20)[["gene_id", "delta_mean_d1_minus_d2", "q_value_bh"]])
    print(top_down.head(20)[["gene_id", "delta_mean_d1_minus_d2", "q_value_bh"]])


if __name__ == "__main__":
    main()
