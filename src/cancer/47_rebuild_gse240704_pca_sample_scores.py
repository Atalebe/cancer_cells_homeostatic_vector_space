#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selected-probe-matrix", required=True)
    ap.add_argument("--outcsv", required=True)
    ap.add_argument("--outjson", required=True)
    ap.add_argument("--n-components", type=int, default=5)
    args = ap.parse_args()

    p = Path(args.selected_probe_matrix)
    if not p.exists():
        raise SystemExit(f"missing input: {p}")

    df = pd.read_parquet(p)

    if "ID_REF" not in df.columns:
        raise SystemExit("selected probe matrix missing ID_REF column")

    sample_cols = [c for c in df.columns if c != "ID_REF"]
    if len(sample_cols) == 0:
        raise SystemExit("no sample columns found in selected probe matrix")

    X = df[sample_cols].to_numpy(dtype=float).T  # samples x probes

    col_means = np.nanmean(X, axis=0)
    inds = np.where(np.isnan(X))
    X[inds] = col_means[inds[1]]

    pca = PCA(n_components=min(args.n_components, X.shape[0], X.shape[1]), svd_solver="auto")
    scores = pca.fit_transform(X)

    out = pd.DataFrame({"sample_id": sample_cols})
    for i in range(scores.shape[1]):
        out[f"PC{i+1}"] = scores[:, i]

    outcsv = Path(args.outcsv)
    outcsv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(outcsv, index=False)

    summary = {
        "selected_probe_matrix": str(p),
        "n_samples": int(X.shape[0]),
        "n_probes": int(X.shape[1]),
        "n_components": int(scores.shape[1]),
        "explained_variance_ratio": {f"PC{i+1}": float(v) for i, v in enumerate(pca.explained_variance_ratio_)},
    }
    with open(args.outjson, "w") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote PCA sample scores:", outcsv)
    print("[ok] wrote summary:", args.outjson)
    print("[info] n_samples:", X.shape[0])
    print("[info] n_probes:", X.shape[1])
    print("[info] explained variance ratios:")
    for i, v in enumerate(pca.explained_variance_ratio_, start=1):
        print(f"  PC{i}: {v:.6f}")


if __name__ == "__main__":
    main()
