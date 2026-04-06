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


PROGRAMS = {
    "antigen_presentation": [
        "B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
        "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5"
    ],
    "ribosomal_translation": [
        "RPL3", "RPL4", "RPL5", "RPL19", "RPLP0",
        "RPS3", "RPS4X", "RPS6", "EEF1A1", "PABPC1", "RACK1"
    ],
    "housekeeping_output": [
        "ACTB", "ACTG1", "EEF1A1", "PABPC1", "HSPA8", "DDX5", "FTL", "RACK1", "GAPDH", "ENO1"
    ],
    "mitochondrial_encoded": [
        "MT-RNR1", "MT-RNR2", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
        "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ATP6", "MT-ATP8"
    ],
    "immune_leaning": [
        "IGHM", "MS4A1", "FCRL3", "FCMR", "BTLA", "THEMIS", "SKAP1", "CRTAM", "KLRB1", "GZMM", "MAL", "FOXO1"
    ],
}


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def strip_version(gene_id: str) -> str:
    return str(gene_id).split(".", 1)[0]


def robust_z(x: np.ndarray) -> np.ndarray:
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    if mad == 0 or np.isnan(mad):
        return np.zeros_like(x, dtype=float)
    return 0.6745 * (x - med) / mad


def mwu_test(a: pd.Series, b: pd.Series) -> dict:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) == 0 or len(b) == 0:
        return {
            "n_a": len(a), "n_b": len(b), "median_a": None, "median_b": None,
            "delta_median_a_minus_b": None, "u_statistic": None, "p_value": None
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
    d2_dir = REPO_ROOT / cfg["results_dir"] / "d2_only_reanalysis"
    out_dir = REPO_ROOT / cfg["results_dir"] / "patient3_d2_1_treatment_sensitivity"
    out_dir.mkdir(parents=True, exist_ok=True)

    state = pd.read_parquet(proc_dir / "state_table.parquet")
    meta = pd.read_parquet(proc_dir / "cell_metadata_registry_aligned.parquet")
    d2 = pd.read_parquet(d2_dir / "d2_state_domains.parquet")[["cell_id", "d2_subdomain"]]

    merged = state.merge(
        meta[["cell_id", "patient_or_donor_id", "treatment_state", "source_class"]],
        on="cell_id", how="left"
    ).merge(d2, on="cell_id", how="left")

    sub = merged.loc[
        (merged["patient_or_donor_id"] == "Patient3") &
        (merged["d2_subdomain"] == "D2_1") &
        (merged["treatment_state"].isin(["treated", "untreated"]))
    ].copy()

    hrsm_rows = []
    for metric in ["H", "S", "M", "R", "phi"]:
        a = sub.loc[sub["treatment_state"] == "treated", metric]
        b = sub.loc[sub["treatment_state"] == "untreated", metric]
        res = mwu_test(a, b)
        res["metric"] = metric
        res["group_a"] = "treated"
        res["group_b"] = "untreated"
        hrsm_rows.append(res)

    hrsm_df = pd.DataFrame(hrsm_rows)[[
        "metric", "group_a", "group_b", "n_a", "n_b",
        "median_a", "median_b", "delta_median_a_minus_b",
        "u_statistic", "p_value"
    ]]
    hrsm_df.to_csv(out_dir / "patient3_d2_1_treatment_hrsm_tests.csv", index=False)

    mat = sparse.load_npz(proc_dir / "counts_sparse.npz").tocsc()
    gene_ids = pd.read_csv(proc_dir / "gene_ids.csv")["gene_id"].astype(str).tolist()
    cell_ids = pd.read_csv(proc_dir / "cell_ids.csv")["cell_id"].astype(str).tolist()

    ann = pd.read_csv(ANNOT_PATH, sep="\t")
    ann["gene_id_stripped"] = ann["gene_id"].astype(str).map(strip_version)

    gene_df = pd.DataFrame({"gene_id": gene_ids})
    gene_df["gene_id_stripped"] = gene_df["gene_id"].astype(str).map(strip_version)
    gene_df = gene_df.merge(
        ann[["gene_id_stripped", "gene_symbol", "is_mito"]],
        on="gene_id_stripped",
        how="left",
    )

    state_map = merged.set_index("cell_id").reindex(cell_ids)
    keep = (
        (state_map["patient_or_donor_id"] == "Patient3") &
        (state_map["d2_subdomain"] == "D2_1") &
        (state_map["treatment_state"].isin(["treated", "untreated"]))
    ).to_numpy()

    mat_sub = mat[:, keep]
    meta_sub = state_map.loc[keep, ["treatment_state"]].reset_index()

    total_counts = np.asarray(mat_sub.sum(axis=0)).ravel().astype(float)
    mito_mask = gene_df["is_mito"].fillna(False).infer_objects(copy=False).to_numpy(dtype=bool)
    mito_counts = np.asarray(mat_sub[mito_mask, :].sum(axis=0)).ravel().astype(float)
    non_mito_counts = total_counts - mito_counts

    score_df = meta_sub.copy()
    score_df["total_counts"] = total_counts
    score_df["mito_counts"] = mito_counts
    score_df["non_mito_counts"] = non_mito_counts
    score_df["mito_fraction"] = mito_counts / np.maximum(total_counts, 1.0)

    for program_name, genes in PROGRAMS.items():
        mask = gene_df["gene_symbol"].isin(genes).to_numpy(dtype=bool)
        if int(mask.sum()) == 0:
            score = np.zeros(len(score_df), dtype=float)
        else:
            raw = np.asarray(mat_sub[mask, :].sum(axis=0)).ravel().astype(float)
            score = robust_z(np.log1p(raw))
        score_df[f"{program_name}_score"] = score

    score_df.to_csv(out_dir / "patient3_d2_1_treatment_program_scores.csv", index=False)

    prog_rows = []
    for metric in [
        "total_counts", "mito_counts", "non_mito_counts", "mito_fraction",
        "antigen_presentation_score", "ribosomal_translation_score",
        "housekeeping_output_score", "mitochondrial_encoded_score", "immune_leaning_score"
    ]:
        a = score_df.loc[score_df["treatment_state"] == "treated", metric]
        b = score_df.loc[score_df["treatment_state"] == "untreated", metric]
        res = mwu_test(a, b)
        res["metric"] = metric
        res["group_a"] = "treated"
        res["group_b"] = "untreated"
        prog_rows.append(res)

    prog_df = pd.DataFrame(prog_rows)[[
        "metric", "group_a", "group_b", "n_a", "n_b",
        "median_a", "median_b", "delta_median_a_minus_b",
        "u_statistic", "p_value"
    ]]
    prog_df.to_csv(out_dir / "patient3_d2_1_treatment_program_overlay_tests.csv", index=False)

    summary = {
        "patient": "Patient3",
        "subset": "D2_1",
        "treated_n": int((sub["treatment_state"] == "treated").sum()),
        "untreated_n": int((sub["treatment_state"] == "untreated").sum()),
    }
    with open(out_dir / "patient3_d2_1_treatment_sensitivity_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))
    print("\n[Patient3 D2_1 HRSM]")
    print(hrsm_df)
    print("\n[Patient3 D2_1 program overlay]")
    print(prog_df)


if __name__ == "__main__":
    main()
