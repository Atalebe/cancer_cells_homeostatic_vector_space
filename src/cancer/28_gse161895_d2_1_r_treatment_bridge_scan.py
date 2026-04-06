#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"
ANNOT_PATH = REPO_ROOT / "data" / "reference" / "annotations" / "ensembl_gene_annotation.tsv"
MIN_DETECTED_CELLS = 10
PROGRESS_EVERY = 500


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def strip_version(gene_id: str) -> str:
    return str(gene_id).split(".", 1)[0]


def main() -> None:
    cfg = read_config()
    proc_dir = REPO_ROOT / cfg["processed_dir"]
    d2_dir = REPO_ROOT / cfg["results_dir"] / "d2_only_reanalysis"
    out_dir = REPO_ROOT / cfg["results_dir"] / "d2_1_r_treatment_bridge"
    out_dir.mkdir(parents=True, exist_ok=True)

    state = pd.read_parquet(proc_dir / "state_table.parquet")[["cell_id", "R"]]
    meta = pd.read_parquet(proc_dir / "cell_metadata_registry_aligned.parquet")[
        ["cell_id", "treatment_state", "patient_or_donor_id", "source_class"]
    ]
    d2 = pd.read_parquet(d2_dir / "d2_state_domains.parquet")[["cell_id", "d2_subdomain"]]

    cells_df = state.merge(meta, on="cell_id", how="left").merge(d2, on="cell_id", how="left")
    cells_df = cells_df.loc[
        (cells_df["d2_subdomain"] == "D2_1") &
        (cells_df["treatment_state"].isin(["treated", "untreated"]))
    ].copy()

    cell_ids = pd.read_csv(proc_dir / "cell_ids.csv")["cell_id"].astype(str)
    gene_ids = pd.read_csv(proc_dir / "gene_ids.csv")["gene_id"].astype(str)
    mat = sparse.load_npz(proc_dir / "counts_sparse.npz").tocsr()

    keep_mask = cell_ids.isin(cells_df["cell_id"]).to_numpy()
    kept_cell_ids = cell_ids[keep_mask].reset_index(drop=True)
    mat_sub = mat[:, keep_mask]

    cell_meta = cells_df.set_index("cell_id").reindex(kept_cell_ids).reset_index()
    treated_mask = (cell_meta["treatment_state"] == "treated").to_numpy()
    untreated_mask = (cell_meta["treatment_state"] == "untreated").to_numpy()
    r_values = cell_meta["R"].to_numpy(dtype=float)

    ann = pd.read_csv(ANNOT_PATH, sep="\t")
    ann["gene_id_stripped"] = ann["gene_id"].astype(str).map(strip_version)

    gene_df = pd.DataFrame({"gene_id": gene_ids})
    gene_df["gene_id_stripped"] = gene_df["gene_id"].astype(str).map(strip_version)
    gene_df = gene_df.merge(
        ann[["gene_id_stripped", "gene_symbol", "gene_biotype", "chromosome"]],
        on="gene_id_stripped",
        how="left",
    )

    # Prefilter genes with too little support
    detected_counts = np.asarray((mat_sub > 0).sum(axis=1)).ravel()
    keep_gene_mask = detected_counts >= MIN_DETECTED_CELLS

    kept_gene_idx = np.where(keep_gene_mask)[0]
    print(f"[info] genes total: {mat_sub.shape[0]}")
    print(f"[info] genes passing MIN_DETECTED_CELLS={MIN_DETECTED_CELLS}: {len(kept_gene_idx)}")

    rows = []
    for j, i in enumerate(kept_gene_idx, start=1):
        x = mat_sub.getrow(i).toarray().ravel().astype(float)
        x_log = np.log1p(x)

        x_treated = x_log[treated_mask]
        x_untreated = x_log[untreated_mask]

        # Treatment contrast
        try:
            u_stat, p_treat = mannwhitneyu(x_treated, x_untreated, alternative="two-sided")
        except ValueError:
            u_stat, p_treat = np.nan, np.nan

        delta = float(np.median(x_treated) - np.median(x_untreated))

        # Skip constant inputs for correlation
        if np.all(x_log == x_log[0]):
            rho, p_r = 0.0, 1.0
        else:
            rho, p_r = spearmanr(x_log, r_values)
            if np.isnan(rho):
                rho, p_r = 0.0, 1.0

        rows.append({
            "gene_id": gene_ids.iloc[i],
            "median_log1p_treated": float(np.median(x_treated)),
            "median_log1p_untreated": float(np.median(x_untreated)),
            "delta_log1p_treated_minus_untreated": delta,
            "u_statistic": float(u_stat) if not np.isnan(u_stat) else np.nan,
            "p_treatment": float(p_treat) if not np.isnan(p_treat) else np.nan,
            "spearman_rho_with_R": float(rho),
            "p_R_correlation": float(p_r),
            "detected_frac_treated": float((x[treated_mask] > 0).mean()),
            "detected_frac_untreated": float((x[untreated_mask] > 0).mean()),
            "n_detected_total": int(detected_counts[i]),
        })

        if j % PROGRESS_EVERY == 0:
            print(f"[info] processed {j} / {len(kept_gene_idx)} kept genes")

    res = pd.DataFrame(rows)
    res["q_treatment_bh"] = multipletests(res["p_treatment"].fillna(1.0), method="fdr_bh")[1]
    res["q_R_correlation_bh"] = multipletests(res["p_R_correlation"].fillna(1.0), method="fdr_bh")[1]

    res = res.merge(gene_df, on="gene_id", how="left")

    res["bridge_direction"] = np.select(
        [
            (res["delta_log1p_treated_minus_untreated"] > 0) & (res["spearman_rho_with_R"] > 0),
            (res["delta_log1p_treated_minus_untreated"] < 0) & (res["spearman_rho_with_R"] < 0),
        ],
        [
            "treated_up_R_positive",
            "treated_down_R_negative",
        ],
        default="discordant_or_null",
    )

    res["bridge_score"] = (
        np.abs(res["delta_log1p_treated_minus_untreated"]) *
        np.abs(res["spearman_rho_with_R"])
    )

    res.to_csv(out_dir / "d2_1_r_treatment_bridge_all_genes.csv", index=False)

    strong = res.loc[
        (res["q_treatment_bh"] < 0.05) &
        (res["q_R_correlation_bh"] < 0.05) &
        (res["bridge_direction"] != "discordant_or_null")
    ].copy()

    strong = strong.sort_values(
        ["bridge_direction", "bridge_score", "q_treatment_bh", "q_R_correlation_bh"],
        ascending=[True, False, True, True],
    )
    strong.to_csv(out_dir / "d2_1_r_treatment_bridge_strong_hits.csv", index=False)

    up = strong.loc[strong["bridge_direction"] == "treated_up_R_positive"].copy()
    down = strong.loc[strong["bridge_direction"] == "treated_down_R_negative"].copy()

    up.to_csv(out_dir / "d2_1_r_treatment_bridge_up_hits.csv", index=False)
    down.to_csv(out_dir / "d2_1_r_treatment_bridge_down_hits.csv", index=False)

    summary = {
        "n_d2_1_cells": int(len(cell_meta)),
        "treated_n": int(treated_mask.sum()),
        "untreated_n": int(untreated_mask.sum()),
        "n_genes_total": int(mat_sub.shape[0]),
        "n_genes_passing_prefilter": int(len(kept_gene_idx)),
        "n_strong_hits": int(len(strong)),
        "n_up_hits": int(len(up)),
        "n_down_hits": int(len(down)),
        "top_up_symbols": up["gene_symbol"].fillna(up["gene_id"]).head(20).tolist(),
        "top_down_symbols": down["gene_symbol"].fillna(down["gene_id"]).head(20).tolist(),
    }

    with open(out_dir / "d2_1_r_treatment_bridge_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))
    print("\n[top up hits]")
    print(up[["gene_id", "gene_symbol", "delta_log1p_treated_minus_untreated", "spearman_rho_with_R", "q_treatment_bh", "q_R_correlation_bh"]].head(20))
    print("\n[top down hits]")
    print(down[["gene_id", "gene_symbol", "delta_log1p_treated_minus_untreated", "spearman_rho_with_R", "q_treatment_bh", "q_R_correlation_bh"]].head(20))


if __name__ == "__main__":
    main()
