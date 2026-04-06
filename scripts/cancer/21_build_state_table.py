import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

from src.common.io import read_table, write_table
from src.common.qc_utils import robust_zscore


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    state_in = processed_dir / "hrsm_proxy_table.csv"

    if not state_in.exists():
        raise FileNotFoundError(f"Missing proxy table: {state_in}")

    df = read_table(state_in)

    df["H"] = robust_zscore(df["H_raw"])
    df["S"] = robust_zscore(df["S_raw"])
    df["M"] = robust_zscore(df["M_raw"])
    df["R"] = robust_zscore(df["R_raw"])

    state_cols = [
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
    ]
    state_df = df[state_cols].copy()

    out_csv = processed_dir / "state_table.csv"
    out_parquet = processed_dir / "state_table.parquet"
    write_table(state_df, out_csv)
    state_df.to_parquet(out_parquet, index=False)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "n_cells": int(state_df.shape[0]),
        "output_csv": str(out_csv),
        "output_parquet": str(out_parquet),
        "group_means": (
            state_df.groupby("population_label")[["H", "S", "M", "R"]]
            .mean()
            .round(4)
            .to_dict(orient="index")
        ),
    }

    with open(processed_dir / "state_table_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
