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


IMMUNE_PROGRAMS = {
    "bcell_lymphoid_like": [
        "IGHM", "MS4A1", "FCRL3", "FCMR", "BTLA", "THEMIS", "SKAP1",
        "CRTAM", "KLRB1", "GZMM", "MAL", "FOXO1", "P2RY10"
    ],
    "antigen_presentation": [
        "B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
        "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5"
    ],
    "stress_immediate_early": [
        "GADD45A", "JUN", "FOS", "DUSP1", "BTG1", "HSPA1A", "HSPA1B"
    ],
    "malignant_output_core": [
        "ACTB", "ACTG1", "EEF1A1", "PABPC1", "HSPA8", "DDX5",
        "RPL3", "RPL4", "RPLP0", "RPS4X", "RPS6", "GAPDH", "ENO1", "FTL"
    ],
}


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


def robust_z(x: np.ndarray) -> np.ndarray:
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    if mad == 0 or np.isnan(mad):
        return np.zeros_like(x, dtype=float)
    return 0.6745 * (x - med) / mad


def main() -> None:
    cfg = read_config()
    proc_dir = REPO_ROOT / cfg["processed_dir"]
    d2_dir = REPO_ROOT / cfg["results_dir"] / "d2_only_reanalysis"
    out_dir = REPO_ROOT / cfg["results_dir"] / "d2_immune_identity_overlay"
    out_dir.mkdir(parents=True, exist_ok=True)

    mat = sparse.load_npz(proc_dir / "counts_sparse.npz").tocsc()
    gene_ids = pd.read_csv(proc_dir / "gene_ids.csv")["gene_id"].astype(str).tolist()
    cell_ids = pd.read_csv(proc_dir / "cell_ids.csv")["cell_id"].astype(str).tolist()

    ann = pd.read_csv(ANNOT_PATH, sep="\t")
    ann["gene_id_stripped"] = ann["gene_id"].astype(str).map(strip_version)

    gene_df = pd.DataFrame({"gene_id": gene_ids})
    gene_df["gene_id_stripped"] = gene_df["gene_id"].astype(str).map(strip_version)
    gene_df = gene_df.merge(
        ann[["gene_id_stripped", "gene_symbol"]],
        on="gene_id_stripped",
        how="left",
    )

    d2 = pd.read_parquet(d2_dir / "d2_state_domains.parquet")[["cell_id", "d2_subdomain"]]
    d2 = d2.set_index("cell_id").reindex(cell_ids).reset_index()
    keep = d2["d2_subdomain"].notna().to_numpy()
    sub = d2.loc[d2["d2_subdomain"].notna(), :].copy()
    mat_d2 = mat[:, keep]

    groups = sorted(sub["d2_subdomain"].astype(str).unique().tolist())
    if len(groups) != 2:
        raise ValueError(f"Expected 2 D2 subdomains, found {groups}")
    group_a, group_b = groups

    score_df = pd.DataFrame({
        "cell_id": sub["cell_id"].tolist(),
        "d2_subdomain": sub["d2_subdomain"].tolist(),
    })

    membership_rows = []
    for program_name, genes in IMMUNE_PROGRAMS.items():
        mask = gene_df["gene_symbol"].isin(genes).to_numpy(dtype=bool)
        present_genes = gene_df.loc[mask, "gene_symbol"].dropna().astype(str).tolist()
        membership_rows.append({
            "program": program_name,
            "n_requested_genes": int(len(genes)),
            "n_present_genes": int(len(present_genes)),
            "present_genes": "; ".join(sorted(set(present_genes))),
        })
        if int(mask.sum()) == 0:
            score = np.zeros(len(score_df), dtype=float)
        else:
            raw = np.asarray(mat_d2[mask, :].sum(axis=0)).ravel().astype(float)
            score = robust_z(np.log1p(raw))
        score_df[f"{program_name}_score"] = score

    score_df.to_csv(out_dir / "d2_immune_identity_cell_scores.csv", index=False)
    pd.DataFrame(membership_rows).to_csv(out_dir / "d2_immune_identity_program_membership.csv", index=False)

    rows = []
    for program_name in IMMUNE_PROGRAMS:
        metric = f"{program_name}_score"
        a = score_df.loc[score_df["d2_subdomain"] == group_a, metric]
        b = score_df.loc[score_df["d2_subdomain"] == group_b, metric]
        res = mwu_test(a, b)
        res["metric"] = metric
        res["group_a"] = group_a
        res["group_b"] = group_b
        rows.append(res)

    summary_df = pd.DataFrame(rows)[[
        "metric", "group_a", "group_b", "n_a", "n_b",
        "median_a", "median_b", "delta_median_a_minus_b",
        "u_statistic", "p_value"
    ]]
    summary_df.to_csv(out_dir / "d2_immune_identity_overlay_tests.csv", index=False)

    with open(out_dir / "d2_immune_identity_overlay_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary_df.to_dict(orient="records"), f, indent=2)

    print(summary_df)


if __name__ == "__main__":
    main()
