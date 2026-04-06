import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

from src.common.io import read_table, write_table


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    tables_dir = Path(cfg["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    state_path = processed_dir / "state_table.csv"
    refined_s_path = tables_dir / "refined_structural_stability_table.csv"
    adjusted_r_path = tables_dir / "adjusted_recoverability_table.csv"

    if not state_path.exists():
        raise FileNotFoundError(f"Missing state table: {state_path}")
    if not refined_s_path.exists():
        raise FileNotFoundError(f"Missing refined structural stability table: {refined_s_path}")
    if not adjusted_r_path.exists():
        raise FileNotFoundError(f"Missing adjusted recoverability table: {adjusted_r_path}")

    state = read_table(state_path)
    sref = read_table(refined_s_path)[["cell_id", "S_refined", "S_refined_raw"]]
    radj = read_table(adjusted_r_path)[["cell_id", "R_adjusted", "R_adjusted_raw"]]

    merged = state.merge(sref, on="cell_id", how="left")
    merged = merged.merge(radj, on="cell_id", how="left")

    merged["S_base"] = merged["S"]
    merged["R_base"] = merged["R"]

    merged["S"] = merged["S_refined"]
    merged["R"] = merged["R_adjusted"]

    out_cols = [
        "cell_id",
        "population_label",
        "group_order",
        "total_counts",
        "detected_genes",
        "zero_fraction",
        "pc1",
        "pc2",
        "pc3",
        "pc4",
        "pc5",
        "H",
        "S",
        "M",
        "R",
        "S_base",
        "R_base",
        "S_refined",
        "R_adjusted",
    ]
    out_df = merged[out_cols].copy()

    out_csv = processed_dir / "state_table_refined.csv"
    out_parquet = processed_dir / "state_table_refined.parquet"
    write_table(out_df, out_csv)
    out_df.to_parquet(out_parquet, index=False)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "n_cells": int(out_df.shape[0]),
        "output_csv": str(out_csv),
        "output_parquet": str(out_parquet),
        "group_means": (
            out_df.groupby("population_label")[["H", "S", "M", "R"]]
            .mean()
            .round(4)
            .to_dict(orient="index")
        ),
        "note": "Promoted refined state table uses S_refined and R_adjusted as default S and R."
    }

    with open(processed_dir / "state_table_refined_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
