#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta-matrix", required=True)
    ap.add_argument("--sample-registry", required=True)
    ap.add_argument("--outparquet", required=True)
    ap.add_argument("--n-probes", type=int, default=5000)
    args = ap.parse_args()

    beta = pd.read_parquet(args.beta_matrix)
    reg = pd.read_csv(args.sample_registry)

    sample_cols = [c for c in beta.columns if c != "ID_REF"]
    if len(sample_cols) == 0:
        raise SystemExit("no sample columns found in beta matrix")

    arr = beta[sample_cols].to_numpy(dtype=np.float32, copy=False)
    variance = np.nanvar(arr, axis=1)
    variance = np.where(np.isfinite(variance), variance, -np.inf)

    order = np.argsort(variance)[::-1]
    keep = order[: args.n_probes]

    out = beta.iloc[keep].copy()
    out.to_parquet(args.outparquet, index=False)

    print(f"[ok] wrote selected probe matrix: {args.outparquet}")
    print(f"[info] selected probes: {len(out)}")
    print(f"[info] sample columns: {len(sample_cols)}")
    print(f"[info] top variance min among kept: {float(np.nanmin(variance[keep]))}")
    print(f"[info] top variance max among kept: {float(np.nanmax(variance[keep]))}")


if __name__ == "__main__":
    main()
