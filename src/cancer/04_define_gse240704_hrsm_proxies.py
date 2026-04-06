from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import NearestNeighbors

from _gse240704_utils import entropy_from_distances, load_config, robust_z


def _parquet_batches(parquet_path: str, batch_size: int = 2048):
    pf = pq.ParquetFile(parquet_path)
    for batch in pf.iter_batches(batch_size=batch_size):
        yield batch.to_pandas()


def _select_top_variable_probes(
    beta_path: str,
    sample_ids: list[str],
    min_non_missing_fraction: float,
    top_n: int,
    batch_size: int = 2048,
) -> pd.DataFrame:
    candidates = []

    batch_no = 0
    for batch_df in _parquet_batches(beta_path, batch_size=batch_size):
        batch_no += 1

        values = batch_df[sample_ids].apply(pd.to_numeric, errors="coerce").astype(np.float32)
        arr = values.to_numpy(dtype=np.float32, copy=False)

        non_missing_fraction = np.isfinite(arr).mean(axis=1)
        variance = np.nanvar(arr, axis=1)

        keep = non_missing_fraction >= min_non_missing_fraction
        if np.any(keep):
            sub = pd.DataFrame(
                {
                    "ID_REF": batch_df.loc[keep, "ID_REF"].astype(str).values,
                    "variance": variance[keep],
                    "non_missing_fraction": non_missing_fraction[keep],
                }
            )
            candidates.append(sub)

        if batch_no % 25 == 0:
            print(f"[info] first pass batch {batch_no}")

    if not candidates:
        raise RuntimeError("No probes survived non-missing filtering in first pass")

    cand = pd.concat(candidates, axis=0, ignore_index=True)
    cand = cand.sort_values("variance", ascending=False).head(top_n).reset_index(drop=True)
    return cand


def _load_selected_probe_matrix(
    beta_path: str,
    sample_ids: list[str],
    selected_ids: set[str],
    batch_size: int = 2048,
) -> pd.DataFrame:
    chunks = []
    batch_no = 0

    for batch_df in _parquet_batches(beta_path, batch_size=batch_size):
        batch_no += 1
        keep = batch_df["ID_REF"].astype(str).isin(selected_ids)
        if keep.any():
            sub = batch_df.loc[keep, ["ID_REF"] + sample_ids].copy()
            chunks.append(sub)

        if batch_no % 25 == 0:
            print(f"[info] second pass batch {batch_no}")

    if not chunks:
        raise RuntimeError("Second pass did not recover any selected probes")

    out = pd.concat(chunks, axis=0, ignore_index=True)
    out = out.drop_duplicates(subset=["ID_REF"]).copy()
    return out


def main(config_path: str) -> None:
    cfg = load_config(config_path)

    beta_path = cfg["outputs"]["beta_matrix"]
    out_path = cfg["outputs"]["hrsm_proxy_table"]

    top_n = int(cfg["qc"]["top_variable_probe_n"])
    n_pcs = int(cfg["qc"]["pca_components"])
    k = int(cfg["qc"]["knn_k"])
    min_non_missing_fraction = float(cfg["qc"]["min_non_missing_fraction_per_probe"])

    header_df = pd.read_parquet(beta_path, columns=None)
    sample_ids = [c for c in header_df.columns if c != "ID_REF"]
    del header_df

    print("[info] first pass, selecting top variable probes")
    top_probe_df = _select_top_variable_probes(
        beta_path=beta_path,
        sample_ids=sample_ids,
        min_non_missing_fraction=min_non_missing_fraction,
        top_n=top_n,
        batch_size=2048,
    )
    selected_ids = set(top_probe_df["ID_REF"].astype(str).tolist())
    print(f"[info] selected probes: {len(selected_ids)}")

    print("[info] second pass, loading selected probe matrix")
    selected = _load_selected_probe_matrix(
        beta_path=beta_path,
        sample_ids=sample_ids,
        selected_ids=selected_ids,
        batch_size=2048,
    )

    selected = selected.set_index("ID_REF")
    selected = selected.loc[selected.index.intersection(top_probe_df["ID_REF"])].copy()

    X = selected.T
    X.index.name = "sample_id"

    X = X.apply(pd.to_numeric, errors="coerce").astype(np.float32)

    imputer = SimpleImputer(strategy="median")
    X_imp = imputer.fit_transform(X)

    max_pcs = min(n_pcs, X_imp.shape[0], X_imp.shape[1])
    pca = PCA(n_components=max_pcs)
    pcs = pca.fit_transform(X_imp)

    reduced_dim = min(5, pcs.shape[1])

    nbrs = NearestNeighbors(n_neighbors=min(k + 1, len(X.index)))
    nbrs.fit(pcs[:, :reduced_dim])
    distances, indices = nbrs.kneighbors(pcs[:, :reduced_dim])

    neighbor_dist = distances[:, 1:]
    mean_neighbor_distance = neighbor_dist.mean(axis=1)
    local_entropy = entropy_from_distances(neighbor_dist)
    s_raw = -pd.Series(local_entropy, index=X.index)

    dmat = pairwise_distances(pcs[:, :reduced_dim])
    np.fill_diagonal(dmat, np.nan)
    residual_distance = np.nanmean(dmat, axis=1)
    r_raw = -pd.Series(residual_distance, index=X.index)

    h_raw = pd.Series(pcs[:, 0], index=X.index)
    if pcs.shape[1] >= 2:
        m_raw = pd.Series(pcs[:, 1], index=X.index)
    else:
        m_raw = pd.Series(np.zeros(len(X.index), dtype=np.float32), index=X.index)

    proxy = pd.DataFrame(
        {
            "sample_id": X.index,
            "H_raw": h_raw.values,
            "S_raw": s_raw.values,
            "M_raw": m_raw.values,
            "R_raw": r_raw.values,
            "mean_neighbor_distance": mean_neighbor_distance,
        }
    )

    proxy["H"] = robust_z(proxy["H_raw"])
    proxy["S"] = robust_z(proxy["S_raw"])
    proxy["M"] = robust_z(proxy["M_raw"])
    proxy["R"] = robust_z(proxy["R_raw"])

    proxy["phi"] = np.sqrt(proxy["H"] ** 2 + proxy["S"] ** 2 + proxy["M"] ** 2 + proxy["R"] ** 2)
    proxy["pca_explained_variance_pc1"] = pca.explained_variance_ratio_[0]
    proxy["pca_explained_variance_pc2"] = (
        pca.explained_variance_ratio_[1] if len(pca.explained_variance_ratio_) > 1 else np.nan
    )
    proxy["selected_probe_count"] = selected.shape[0]

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    proxy.to_csv(out_path, index=False)

    print(f"[ok] wrote HRSM proxy table: {out_path}")
    print(f"[info] selected probe count: {selected.shape[0]}")
    print("[info] explained variance ratios:")
    for i, v in enumerate(pca.explained_variance_ratio_[:5], start=1):
        print(f"  PC{i}: {v:.6f}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    args = parser.parse_args()
    main(args.config)
