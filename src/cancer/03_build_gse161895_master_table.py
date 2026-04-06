#!/usr/bin/env python3

from __future__ import annotations

import gzip
import json
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import sparse


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"

SPECIAL_HTSEQ_ROWS = {
    "__no_feature",
    "__ambiguous",
    "__too_low_aQual",
    "__not_aligned",
    "__alignment_not_unique",
}


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def extract_tar_if_needed(raw_dir: Path) -> list[Path]:
    tar_path = raw_dir / "GSE161895_RAW.tar"
    extract_dir = raw_dir / "GSE161895_raw_extracted"
    extract_dir.mkdir(parents=True, exist_ok=True)

    existing = sorted(extract_dir.glob("*.htseq-count.txt.gz"))
    if existing:
        return existing

    if not tar_path.exists():
        raise FileNotFoundError(f"Missing tar archive: {tar_path}")

    print(f"[info] extracting {tar_path} to {extract_dir}")
    with tarfile.open(tar_path, "r") as tar:
        tar.extractall(path=extract_dir)

    files = sorted(extract_dir.glob("*.htseq-count.txt.gz"))
    if not files:
        raise RuntimeError("No htseq-count files found after extraction.")
    return files


def build_gene_index(first_file: Path) -> tuple[list[str], dict[str, int]]:
    with gzip.open(first_file, "rt", encoding="utf-8", errors="ignore") as f:
        df = pd.read_csv(f, sep="\t", header=None, names=["gene_id", "count"])
    df = df[~df["gene_id"].isin(SPECIAL_HTSEQ_ROWS)].copy()
    genes = df["gene_id"].astype(str).tolist()
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    return genes, gene_to_idx


def parse_one_count_file_sparse(path: Path, gene_to_idx: dict[str, int]) -> tuple[np.ndarray, np.ndarray]:
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        df = pd.read_csv(f, sep="\t", header=None, names=["gene_id", "count"])

    df = df[~df["gene_id"].isin(SPECIAL_HTSEQ_ROWS)].copy()
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0).astype(np.int32)
    df = df[df["count"] > 0].copy()

    rows = df["gene_id"].map(gene_to_idx).astype(np.int32).values
    vals = df["count"].values.astype(np.int32)
    return rows, vals


def build_sparse_matrix(files: list[Path], gene_to_idx: dict[str, int]) -> tuple[sparse.csc_matrix, list[str]]:
    n_genes = len(gene_to_idx)
    n_cells = len(files)

    row_idx_all = []
    col_idx_all = []
    data_all = []
    cell_ids = []

    for j, path in enumerate(files):
        rows, vals = parse_one_count_file_sparse(path, gene_to_idx)
        cols = np.full(shape=len(rows), fill_value=j, dtype=np.int32)

        row_idx_all.append(rows)
        col_idx_all.append(cols)
        data_all.append(vals)

        cell_id = path.name.replace(".htseq-count.txt.gz", "")
        cell_ids.append(cell_id)

        if (j + 1) % 250 == 0:
            print(f"[info] parsed {j + 1} files")

    row_idx = np.concatenate(row_idx_all)
    col_idx = np.concatenate(col_idx_all)
    data = np.concatenate(data_all)

    mat = sparse.coo_matrix((data, (row_idx, col_idx)), shape=(n_genes, n_cells), dtype=np.int32).tocsc()
    return mat, cell_ids


def build_cell_name_map(cell_ids: list[str]) -> pd.DataFrame:
    rows = []
    for cell_id in cell_ids:
        gsm = cell_id.split("_", 1)[0] if "_" in cell_id else cell_id
        rows.append({"cell_id": cell_id, "gsm": gsm})
    return pd.DataFrame(rows)


def main() -> None:
    cfg = read_config()
    raw_dir = REPO_ROOT / cfg["raw_dir"]
    meta_dir = REPO_ROOT / cfg["metadata_dir"]
    proc_dir = REPO_ROOT / cfg["processed_dir"]
    proc_dir.mkdir(parents=True, exist_ok=True)

    files = extract_tar_if_needed(raw_dir)
    print(f"[info] found {len(files)} htseq-count files")

    genes, gene_to_idx = build_gene_index(files[0])
    mat, cell_ids = build_sparse_matrix(files, gene_to_idx)

    reg = pd.read_parquet(meta_dir / "cell_metadata_registry.parquet").copy()
    cell_map = build_cell_name_map(cell_ids)

    if "gsm" in reg.columns:
        reg["gsm"] = reg["gsm"].astype(str)
    cell_map["gsm"] = cell_map["gsm"].astype(str)

    merged_meta = cell_map.merge(reg, on="gsm", how="left", suffixes=("", "_registry"))

    sparse.save_npz(proc_dir / "counts_sparse.npz", mat)
    pd.DataFrame({"gene_id": genes}).to_csv(proc_dir / "gene_ids.csv", index=False)
    pd.DataFrame({"cell_id": cell_ids}).to_csv(proc_dir / "cell_ids.csv", index=False)
    merged_meta.to_parquet(proc_dir / "cell_metadata_registry_aligned.parquet")
    merged_meta.to_csv(proc_dir / "cell_metadata_registry_aligned.csv", index=False)

    summary = {
        "n_genes": int(mat.shape[0]),
        "n_cells_expression": int(mat.shape[1]),
        "n_nonzero_entries": int(mat.nnz),
        "sparsity": float(1.0 - (mat.nnz / (mat.shape[0] * mat.shape[1]))),
        "n_htseq_files": int(len(files)),
        "n_cells_with_registry_match": int(merged_meta.drop(columns=["cell_id", "gsm"]).notna().any(axis=1).sum()),
        "n_cells_without_registry_match": int((~merged_meta.drop(columns=["cell_id", "gsm"]).notna().any(axis=1)).sum()),
        "n_registry_rows": int(len(reg)),
        "n_unique_gsms_expression": int(cell_map["gsm"].nunique()),
        "n_unique_gsms_registry": int(reg["gsm"].nunique()) if "gsm" in reg.columns else None,
    }

    with open(proc_dir / "master_table_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote sparse realization:")
    print(proc_dir / "counts_sparse.npz")
    print(proc_dir / "gene_ids.csv")
    print(proc_dir / "cell_ids.csv")
    print(proc_dir / "cell_metadata_registry_aligned.parquet")
    print(proc_dir / "master_table_summary.json")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
