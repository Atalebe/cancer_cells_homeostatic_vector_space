#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from scipy.stats import linregress


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def summarize_variance_scaling(mat_sub, label: str) -> dict:
    dense = mat_sub.toarray().astype(float)
    means = dense.mean(axis=1)
    vars_ = dense.var(axis=1, ddof=1)

    mask = (means > 0) & (vars_ > 0)
    log_mean = np.log10(means[mask])
    log_var = np.log10(vars_[mask])

    fit = linregress(log_mean, log_var)

    return {
        "group": label,
        "n_genes_used": int(mask.sum()),
        "slope_logvar_vs_logmean": float(fit.slope),
        "intercept": float(fit.intercept),
        "r_value": float(fit.rvalue),
        "r_squared": float(fit.rvalue ** 2),
        "p_value": float(fit.pvalue),
        "stderr": float(fit.stderr),
        "mean_log_mean": float(log_mean.mean()),
        "mean_log_var": float(log_var.mean()),
    }, pd.DataFrame({
        "mean_expr": means,
        "variance_expr": vars_,
        "used_in_fit": mask,
    })


def main() -> None:
    cfg = read_config()
    proc_dir = REPO_ROOT / cfg["processed_dir"]
    d2_dir = REPO_ROOT / cfg["results_dir"] / "d2_only_reanalysis"
    out_dir = REPO_ROOT / cfg["results_dir"] / "d2_1_variance_scaling"
    out_dir.mkdir(parents=True, exist_ok=True)

    state = pd.read_parquet(proc_dir / "state_table.parquet")[["cell_id"]]
    meta = pd.read_parquet(proc_dir / "cell_metadata_registry_aligned.parquet")[["cell_id", "treatment_state"]]
    d2 = pd.read_parquet(d2_dir / "d2_state_domains.parquet")[["cell_id", "d2_subdomain"]]
    cells_df = state.merge(meta, on="cell_id", how="left").merge(d2, on="cell_id", how="left")

    cell_ids = pd.read_csv(proc_dir / "cell_ids.csv")["cell_id"].astype(str)
    mat = sparse.load_npz(proc_dir / "counts_sparse.npz").tocsc()

    cell_meta = cells_df.set_index("cell_id").reindex(cell_ids).reset_index()

    masks = {
        "D2_1_treated": ((cell_meta["d2_subdomain"] == "D2_1") & (cell_meta["treatment_state"] == "treated")).to_numpy(),
        "D2_1_untreated": ((cell_meta["d2_subdomain"] == "D2_1") & (cell_meta["treatment_state"] == "untreated")).to_numpy(),
        "D2_2_all": (cell_meta["d2_subdomain"] == "D2_2").to_numpy(),
    }

    summary_rows = []
    for label, mask in masks.items():
        if mask.sum() < 10:
            continue
        group_summary, group_df = summarize_variance_scaling(mat[:, mask], label)
        summary_rows.append(group_summary)
        group_df.to_csv(out_dir / f"{label}_gene_mean_variance.csv", index=False)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_dir / "variance_scaling_summary.csv", index=False)

    with open(out_dir / "variance_scaling_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary_rows, f, indent=2)

    print(summary)


if __name__ == "__main__":
    main()
