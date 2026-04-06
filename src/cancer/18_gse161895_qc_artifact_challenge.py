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
    out_dir = REPO_ROOT / cfg["results_dir"] / "qc_artifact_challenge"
    out_dir.mkdir(parents=True, exist_ok=True)

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

    mito_mask = gene_df["is_mito"].fillna(False).infer_objects(copy=False).to_numpy(dtype=bool)
    ribo_mask = gene_df["gene_symbol"].fillna("").astype(str).str.match(r"^(RPL|RPS|MRPL|MRPS)").to_numpy(dtype=bool)

    total_counts = np.asarray(mat.sum(axis=0)).ravel().astype(float)
    detected_genes = np.asarray((mat > 0).sum(axis=0)).ravel().astype(float)
    mito_counts = np.asarray(mat[mito_mask, :].sum(axis=0)).ravel().astype(float)
    ribo_counts = np.asarray(mat[ribo_mask, :].sum(axis=0)).ravel().astype(float)

    mito_fraction = mito_counts / np.maximum(total_counts, 1.0)
    ribo_fraction = ribo_counts / np.maximum(total_counts, 1.0)

    qc_df = pd.DataFrame({
        "cell_id": cell_ids,
        "total_counts": total_counts,
        "detected_genes": detected_genes,
        "mito_counts": mito_counts,
        "mito_fraction": mito_fraction,
        "ribo_counts": ribo_counts,
        "ribo_fraction": ribo_fraction,
        "log10_total_counts": np.log10(np.maximum(total_counts, 1.0)),
        "log10_detected_genes": np.log10(np.maximum(detected_genes, 1.0)),
    })

    domains = pd.read_parquet(state_dir / "state_domains.parquet")[["cell_id", "state_domain"]]
    qc_df = qc_df.merge(domains, on="cell_id", how="left")
    qc_df.to_csv(out_dir / "cell_qc_metrics.csv", index=False)

    rows = []
    metrics = [
        "total_counts",
        "detected_genes",
        "mito_counts",
        "mito_fraction",
        "ribo_counts",
        "ribo_fraction",
        "log10_total_counts",
        "log10_detected_genes",
    ]

    for metric in metrics:
        a = qc_df.loc[qc_df["state_domain"] == "D1", metric]
        b = qc_df.loc[qc_df["state_domain"] == "D2", metric]
        res = mwu_test(a, b)
        res["metric"] = metric
        rows.append(res)

    summary_df = pd.DataFrame(rows)[[
        "metric", "n_a", "n_b", "median_a", "median_b",
        "delta_median_a_minus_b", "u_statistic", "p_value"
    ]]
    summary_df.to_csv(out_dir / "d1_vs_d2_qc_tests.csv", index=False)

    # challenge subset, exclude the lowest-depth tail globally
    total_cut = float(np.quantile(qc_df["total_counts"], 0.10))
    gene_cut = float(np.quantile(qc_df["detected_genes"], 0.10))
    mito_cut = float(np.quantile(qc_df["mito_fraction"], 0.90))

    challenge_subset = qc_df[
        (qc_df["total_counts"] >= total_cut) &
        (qc_df["detected_genes"] >= gene_cut) &
        (qc_df["mito_fraction"] <= mito_cut)
    ].copy()

    challenge_subset.to_csv(out_dir / "challenge_subset_cells.csv", index=False)

    challenge_rows = []
    for metric in metrics:
        a = challenge_subset.loc[challenge_subset["state_domain"] == "D1", metric]
        b = challenge_subset.loc[challenge_subset["state_domain"] == "D2", metric]
        res = mwu_test(a, b)
        res["metric"] = metric
        challenge_rows.append(res)

    challenge_df = pd.DataFrame(challenge_rows)[[
        "metric", "n_a", "n_b", "median_a", "median_b",
        "delta_median_a_minus_b", "u_statistic", "p_value"
    ]]
    challenge_df.to_csv(out_dir / "d1_vs_d2_qc_tests_challenge_subset.csv", index=False)

    summary = {
        "n_cells_total": int(len(qc_df)),
        "n_cells_challenge_subset": int(len(challenge_subset)),
        "challenge_thresholds": {
            "total_counts_ge": total_cut,
            "detected_genes_ge": gene_cut,
            "mito_fraction_le": mito_cut,
        },
    }
    with open(out_dir / "qc_artifact_challenge_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(summary_df)
    print("\n[challenge subset]")
    print(challenge_df)


if __name__ == "__main__":
    main()
