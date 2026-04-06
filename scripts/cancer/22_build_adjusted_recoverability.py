import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src.common.io import read_table, write_table
from src.common.qc_utils import robust_zscore


def fit_linear_residual(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    X_design = np.column_stack([np.ones(len(X)), X])
    beta, *_ = np.linalg.lstsq(X_design, y, rcond=None)
    fitted = X_design @ beta
    resid = y - fitted
    return resid


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--use-full-metabolic-panel", action="store_true")
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    tables_dir = Path(cfg["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    state_path = processed_dir / "state_table.csv"
    metabolic_path = tables_dir / "metabolic_overlay_table.csv"
    refined_s_path = tables_dir / "refined_structural_stability_table.csv"

    if not state_path.exists():
        raise FileNotFoundError(f"Missing state table: {state_path}")
    if not metabolic_path.exists():
        raise FileNotFoundError(f"Missing metabolic overlay table: {metabolic_path}")

    state = read_table(state_path)
    metabolic = read_table(metabolic_path)

    keep_cols = ["cell_id", "glycolysis_score", "oxphos_score", "metabolic_balance_score"]
    merged = state.merge(metabolic[keep_cols], on="cell_id", how="left")

    if args.use_full_metabolic_panel:
        X_cols = ["glycolysis_score", "oxphos_score", "metabolic_balance_score"]
    else:
        X_cols = ["metabolic_balance_score"]

    y = merged["R"].values.astype(float)
    X = merged[X_cols].values.astype(float)

    resid = fit_linear_residual(y, X)
    merged["R_adjusted_raw"] = resid
    merged["R_adjusted"] = robust_zscore(pd.Series(resid))

    merged["R_minus_R_adjusted"] = merged["R"] - merged["R_adjusted"]

    if refined_s_path.exists():
        refined_s = read_table(refined_s_path)[["cell_id", "S_refined", "S_refined_raw"]]
        merged = merged.merge(refined_s, on="cell_id", how="left")
    else:
        merged["S_refined"] = np.nan
        merged["S_refined_raw"] = np.nan

    out_path = tables_dir / "adjusted_recoverability_table.csv"
    write_table(merged, out_path)

    group_summary = (
        merged.groupby("population_label")[
            ["R", "R_adjusted", "metabolic_balance_score", "glycolysis_score", "oxphos_score", "H", "M"]
        ]
        .agg(["mean", "median", "std"])
        .round(4)
    )
    group_summary.columns = ["_".join(col).strip("_") for col in group_summary.columns]
    group_summary = group_summary.reset_index()

    group_summary_path = tables_dir / "adjusted_recoverability_group_summary.csv"
    write_table(group_summary, group_summary_path)

    corr_cols = [
        "H", "S", "M", "R", "R_adjusted",
        "glycolysis_score", "oxphos_score", "metabolic_balance_score"
    ]
    if "S_refined" in merged.columns:
        corr_cols.append("S_refined")

    corr = merged[corr_cols].corr().round(4).reset_index()
    corr_path = tables_dir / "adjusted_recoverability_correlations.csv"
    write_table(corr, corr_path)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "regressors_used": X_cols,
        "outputs": {
            "table": str(out_path),
            "group_summary": str(group_summary_path),
            "correlations": str(corr_path),
        },
        "group_means": (
            merged.groupby("population_label")[["R", "R_adjusted"]]
            .mean()
            .round(4)
            .to_dict(orient="index")
        ),
    }

    with open(tables_dir / "adjusted_recoverability_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
