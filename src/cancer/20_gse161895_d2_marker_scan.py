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
    d2_dir = REPO_ROOT / cfg["results_dir"] / "d2_only_reanalysis"
    out_dir = REPO_ROOT / cfg["results_dir"] / "d2_marker_scan"
    out_dir.mkdir(parents=True, exist_ok=True)

    mat = sparse.load_npz(proc_dir / "counts_sparse.npz").tocsc()
    gene_ids = pd.read_csv(proc_dir / "gene_ids.csv")["gene_id"].astype(str).tolist()
    cell_ids = pd.read_csv(proc_dir / "cell_ids.csv")["cell_id"].astype(str).tolist()

    d2 = pd.read_parquet(d2_dir / "d2_state_domains.parquet")[["cell_id", "d2_subdomain"]]
    d2 = d2.set_index("cell_id").reindex(cell_ids)

    keep_cells = d2["d2_subdomain"].notna().to_numpy()
    sub_labels = d2.loc[d2["d2_subdomain"].notna(), "d2_subdomain"].astype(str).tolist()
    mat_d2 = mat[:, keep_cells]

    groups = sorted(pd.Series(sub_labels).unique().tolist())
    if len(groups) != 2:
        raise ValueError(f"Expected exactly 2 D2 subdomains, found {groups}")

    group_a, group_b = groups
    idx_a = np.where(np.array(sub_labels) == group_a)[0]
    idx_b = np.where(np.array(sub_labels) == group_b)[0]

    mat_a = mat_d2[:, idx_a]
    mat_b = mat_d2[:, idx_b]

    mean_a = np.asarray(mat_a.mean(axis=1)).ravel()
    mean_b = np.asarray(mat_b.mean(axis=1)).ravel()
    det_a = np.asarray((mat_a > 0).mean(axis=1)).ravel()
    det_b = np.asarray((mat_b > 0).mean(axis=1)).ravel()

    keep = (det_a >= 0.05) | (det_b >= 0.05)
    keep_idx = np.where(keep)[0]

    rows = []
    pvals = []

    for i in keep_idx:
        a = mat_a[i, :].toarray().ravel()
        b = mat_b[i, :].toarray().ravel()

        try:
            u, p = mannwhitneyu(a, b, alternative="two-sided")
        except Exception:
            u, p = np.nan, 1.0

        rows.append({
            "gene_id": gene_ids[i],
            "mean_a": float(mean_a[i]),
            "mean_b": float(mean_b[i]),
            "delta_mean_a_minus_b": float(mean_a[i] - mean_b[i]),
            "detected_frac_a": float(det_a[i]),
            "detected_frac_b": float(det_b[i]),
            "delta_detected_frac_a_minus_b": float(det_a[i] - det_b[i]),
            "u_statistic": float(u) if not np.isnan(u) else None,
            "p_value": float(p),
            "group_a": group_a,
            "group_b": group_b,
        })
        pvals.append(p)

    out = pd.DataFrame(rows)
    out["q_value_bh"] = bh_adjust(np.array(pvals))
    out["abs_delta_mean"] = out["delta_mean_a_minus_b"].abs()
    out["is_ercc"] = out["gene_id"].astype(str).str.startswith("ERCC-")

    out = out.sort_values(["q_value_bh", "abs_delta_mean"], ascending=[True, False]).reset_index(drop=True)
    out.to_csv(out_dir / "d2_subdomain_marker_scan_full.csv", index=False)

    clean = out.loc[~out["is_ercc"]].copy()
    clean.to_csv(out_dir / "d2_subdomain_marker_scan_full_clean.csv", index=False)

    high = clean.sort_values(["delta_mean_a_minus_b", "q_value_bh"], ascending=[False, True]).head(200)
    low = clean.sort_values(["delta_mean_a_minus_b", "q_value_bh"], ascending=[True, True]).head(200)

    high.to_csv(out_dir / "d2_1_high_top200_clean.csv", index=False)
    low.to_csv(out_dir / "d2_1_low_top200_clean.csv", index=False)

    summary = {
        "group_a": group_a,
        "group_b": group_b,
        "n_genes_tested": int(len(out)),
        "n_genes_clean": int(len(clean)),
        "n_group_a_cells": int(len(idx_a)),
        "n_group_b_cells": int(len(idx_b)),
    }
    with open(out_dir / "d2_subdomain_marker_scan_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))
    print(high.head(20)[["gene_id", "delta_mean_a_minus_b", "q_value_bh"]])
    print(low.head(20)[["gene_id", "delta_mean_a_minus_b", "q_value_bh"]])


if __name__ == "__main__":
    main()
