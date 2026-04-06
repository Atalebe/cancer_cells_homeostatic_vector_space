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

ANNOT_PATH = REPO_ROOT / "data" / "reference" / "annotations" / "ensembl_gene_annotation.tsv"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def strip_version(gene_id: str) -> str:
    return str(gene_id).split(".", 1)[0]


def mwu_test(a: pd.Series, b: pd.Series) -> dict:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) == 0 or len(b) == 0:
        return {
            "n_a": len(a),
            "n_b": len(b),
            "median_a": None,
            "median_b": None,
            "delta_median_a_minus_b": None,
            "u_statistic": None,
            "p_value": None,
        }
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(a.median()),
        "median_b": float(b.median()),
        "delta_median_a_minus_b": float(a.median() - b.median()),
        "u_statistic": float(u),
        "p_value": float(p),
    }


def main() -> None:
    cfg = read_config()
    proc_dir = REPO_ROOT / cfg["processed_dir"]
    state_dir = REPO_ROOT / cfg["results_dir"] / "state_domains"
    out_dir = REPO_ROOT / cfg["results_dir"] / "mito_fraction_followup"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not ANNOT_PATH.exists():
        raise FileNotFoundError(f"Missing annotation TSV: {ANNOT_PATH}")

    mat = sparse.load_npz(proc_dir / "counts_sparse.npz").tocsc()
    gene_ids = pd.read_csv(proc_dir / "gene_ids.csv")["gene_id"].astype(str).tolist()
    cell_ids = pd.read_csv(proc_dir / "cell_ids.csv")["cell_id"].astype(str).tolist()

    ann = pd.read_csv(ANNOT_PATH, sep="\t")
    ann["gene_id_stripped"] = ann["gene_id"].astype(str).map(strip_version)

    gene_df = pd.DataFrame({"gene_id": gene_ids})
    gene_df["gene_id_stripped"] = gene_df["gene_id"].astype(str).map(strip_version)
    gene_df = gene_df.merge(
        ann[["gene_id_stripped", "gene_symbol", "gene_biotype", "chromosome", "is_mito"]],
        on="gene_id_stripped",
        how="left",
    )

    mito_mask = gene_df["is_mito"].fillna(False).to_numpy(dtype=bool)

    total_counts = np.asarray(mat.sum(axis=0)).ravel().astype(float)
    mito_counts = np.asarray(mat[mito_mask, :].sum(axis=0)).ravel().astype(float)
    mito_frac = mito_counts / np.maximum(total_counts, 1.0)

    mito_df = pd.DataFrame({
        "cell_id": cell_ids,
        "total_counts": total_counts,
        "mito_counts": mito_counts,
        "mito_fraction": mito_frac,
    })

    domains = pd.read_parquet(state_dir / "state_domains.parquet")[["cell_id", "state_domain"]]
    mito_df = mito_df.merge(domains, on="cell_id", how="left")

    mito_df.to_csv(out_dir / "cell_mito_fraction.csv", index=False)

    a = mito_df.loc[mito_df["state_domain"] == "D1", "mito_fraction"]
    b = mito_df.loc[mito_df["state_domain"] == "D2", "mito_fraction"]
    test = mwu_test(a, b)
    test["metric"] = "mito_fraction"
    summary_df = pd.DataFrame([test])
    summary_df.to_csv(out_dir / "d1_vs_d2_mito_fraction_test.csv", index=False)

    summary = {
        "n_mito_genes_annotated": int(mito_mask.sum()),
        "n_cells": int(len(mito_df)),
        "d1_vs_d2_test": test,
    }
    with open(out_dir / "mito_fraction_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
