#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def robust_z(x: pd.Series) -> pd.Series:
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    if mad == 0 or np.isnan(mad):
        return pd.Series(np.zeros(len(x)), index=x.index)
    return 0.6745 * (x - med) / mad


def safe_match_normal_reference(meta: pd.DataFrame, expr_cols: list[str]) -> list[str]:
    candidate_cols = [c for c in meta.columns if c not in {"cell_id", "gsm"}]
    ref_cells = []

    for c in candidate_cols:
        vals = meta[c].astype(str)
        mask = vals.str.contains("normal", case=False, na=False) | vals.str.contains("donor", case=False, na=False)
        if mask.any():
            ref_cells.extend(meta.loc[mask, "cell_id"].astype(str).tolist())

    ref_cells = sorted(set([x for x in ref_cells if x in expr_cols]))
    return ref_cells


def main() -> None:
    cfg = read_config()
    proc_dir = REPO_ROOT / cfg["processed_dir"]
    out_dir = REPO_ROOT / cfg["results_dir"] / "proxies"
    out_dir.mkdir(parents=True, exist_ok=True)

    mat = sparse.load_npz(proc_dir / "counts_sparse.npz").tocsc()
    gene_ids = pd.read_csv(proc_dir / "gene_ids.csv")["gene_id"].astype(str).tolist()
    cell_ids = pd.read_csv(proc_dir / "cell_ids.csv")["cell_id"].astype(str).tolist()
    meta = pd.read_parquet(proc_dir / "cell_metadata_registry_aligned.parquet")

    detected_per_gene = np.asarray((mat > 0).sum(axis=1)).ravel()
    keep_genes = detected_per_gene >= cfg["filters"]["min_cells_per_gene"]
    mat_f = mat[keep_genes, :]

    # library-size normalize and log1p on nonzero entries only
    total_counts = np.asarray(mat.sum(axis=0)).ravel().astype(np.float64)
    detected_genes = np.asarray((mat > 0).sum(axis=0)).ravel().astype(np.int32)

    scale_factor = cfg["normalization"]["log1p_scale_factor"]
    col_scale = scale_factor / np.maximum(total_counts, 1.0)

    norm_mat = mat_f.astype(np.float32).tocsc(copy=True)
    for j in range(norm_mat.shape[1]):
        start, end = norm_mat.indptr[j], norm_mat.indptr[j + 1]
        norm_mat.data[start:end] *= col_scale[j]
    norm_mat.data = np.log1p(norm_mat.data)

    n_components = min(cfg["normalization"]["n_pcs"], norm_mat.shape[0] - 1, norm_mat.shape[1] - 1, 50)
    n_components = max(n_components, 2)

    svd = TruncatedSVD(n_components=n_components, random_state=cfg["seed"])
    pcs = svd.fit_transform(norm_mat.T)

    H = robust_z(pd.Series(np.log1p(total_counts), index=cell_ids)) + robust_z(pd.Series(detected_genes, index=cell_ids))

    k = min(cfg["knn"]["k_stability_default"] + 1, len(cell_ids))
    nbrs = NearestNeighbors(n_neighbors=k)
    nbrs.fit(pcs[:, : min(10, pcs.shape[1])])
    dists, _ = nbrs.kneighbors(pcs[:, : min(10, pcs.shape[1])])
    mean_knn_dist = pd.Series(dists[:, 1:].mean(axis=1), index=cell_ids)
    S = -robust_z(mean_knn_dist)

    M = robust_z(pd.Series(pcs[:, 0], index=cell_ids))

    ref_cells = safe_match_normal_reference(meta, cell_ids)
    if not ref_cells:
        ref_cells = cell_ids[: max(10, min(50, len(cell_ids) // 20))]

    pc_df = pd.DataFrame(pcs[:, : min(10, pcs.shape[1])], index=cell_ids)
    ref_centroid = pc_df.loc[ref_cells].mean(axis=0)
    dist_to_ref = ((pc_df - ref_centroid) ** 2).sum(axis=1).pow(0.5)
    R = -robust_z(dist_to_ref)

    proxy_df = pd.DataFrame({
        "cell_id": cell_ids,
        "H": H.values,
        "S": S.values,
        "M": M.values,
        "R": R.values,
        "log_total_counts": np.log1p(total_counts),
        "detected_genes": detected_genes,
    })

    proxy_df.to_csv(proc_dir / "hrsm_proxy_table.csv", index=False)

    summary = {
        "n_genes_total": int(mat.shape[0]),
        "n_genes_after_filter": int(mat_f.shape[0]),
        "n_cells": int(mat.shape[1]),
        "explained_variance_ratio_first5": svd.explained_variance_ratio_[:5].tolist(),
        "n_reference_cells_for_R": int(len(ref_cells)),
    }
    with open(out_dir / "proxy_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote hrsm_proxy_table.csv and proxy summary")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
