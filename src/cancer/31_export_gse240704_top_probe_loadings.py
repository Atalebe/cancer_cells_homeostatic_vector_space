#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() == ".parquet":
        return pd.read_parquet(p)
    if p.suffix.lower() == ".csv":
        return pd.read_csv(p)
    raise SystemExit(f"unsupported input format: {path}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selected-probe-matrix", required=True)
    ap.add_argument("--sample-annotations-curated", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--top-n", type=int, default=100)
    ap.add_argument("--n-components", type=int, default=5)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    mat = read_table(args.selected_probe_matrix)
    ann = read_table(args.sample_annotations_curated)

    if "ID_REF" not in mat.columns:
        raise SystemExit("selected probe matrix must contain ID_REF")

    sample_cols = [c for c in mat.columns if c != "ID_REF"]
    X = mat[sample_cols].to_numpy(dtype=np.float32).T

    if np.isnan(X).any():
        col_means = np.nanmean(X, axis=0)
        inds = np.where(np.isnan(X))
        X[inds] = col_means[inds[1]]

    n_components = min(args.n_components, X.shape[0], X.shape[1])
    pca = PCA(n_components=n_components, svd_solver="full", random_state=0)
    pca.fit(X)

    probe_ids = mat["ID_REF"].astype(str).to_numpy()
    rows = []
    for i in range(n_components):
        load = pca.components_[i]
        abs_order = np.argsort(np.abs(load))[::-1][: args.top_n]
        for rank, idx in enumerate(abs_order, start=1):
            rows.append(
                {
                    "pc": f"PC{i+1}",
                    "rank_by_abs_loading": rank,
                    "ID_REF": probe_ids[idx],
                    "loading": float(load[idx]),
                    "abs_loading": float(abs(load[idx])),
                }
            )

    out = pd.DataFrame(rows)
    out.to_csv(outdir / "top_probe_loadings_by_pc.csv", index=False)

    ev = pd.DataFrame(
        {
            "pc": [f"PC{i+1}" for i in range(n_components)],
            "explained_variance_ratio": pca.explained_variance_ratio_,
        }
    )
    ev.to_csv(outdir / "pca_explained_variance_rebuilt.csv", index=False)

    print(f"[ok] wrote top probe loadings: {outdir / 'top_probe_loadings_by_pc.csv'}")
    print(f"[ok] wrote explained variance table: {outdir / 'pca_explained_variance_rebuilt.csv'}")
    print(f"[info] n_components: {n_components}")
    print(f"[info] n_probes: {len(probe_ids)}")
    print(f"[info] n_samples: {len(sample_cols)}")


if __name__ == "__main__":
    main()
