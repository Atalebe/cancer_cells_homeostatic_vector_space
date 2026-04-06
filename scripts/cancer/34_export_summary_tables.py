import argparse
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
    if not state_path.exists():
        raise FileNotFoundError(f"Missing state table: {state_path}")

    df = read_table(state_path)

    group_summary = (
        df.groupby("population_label")[["H", "S", "M", "R", "total_counts", "detected_genes", "zero_fraction"]]
        .agg(["mean", "median", "std"])
        .round(4)
    )
    group_summary.columns = ["_".join(col).strip("_") for col in group_summary.columns]
    group_summary = group_summary.reset_index()

    pairwise_summary = df[["H", "S", "M", "R"]].corr().round(4).reset_index()

    write_table(group_summary, tables_dir / "group_state_summary.csv")
    write_table(pairwise_summary, tables_dir / "state_correlation_matrix.csv")

    print(f"[ok] wrote {tables_dir / 'group_state_summary.csv'}")
    print(f"[ok] wrote {tables_dir / 'state_correlation_matrix.csv'}")


if __name__ == "__main__":
    main()
