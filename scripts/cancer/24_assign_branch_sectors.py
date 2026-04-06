import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src.common.io import read_table, write_table


def assign_sector(row, med):
    high_H = row["H"] >= med["H"]
    high_M = row["M"] >= med["M"]
    high_S = row["S"] >= med["S"]
    high_R = row["R"] >= med["R"]

    low_H = row["H"] < med["H"]
    low_M = row["M"] < med["M"]
    low_S = row["S"] < med["S"]
    low_R = row["R"] < med["R"]

    if high_H and high_M and low_R and low_S:
        return "unstable_committed_malignant"

    if high_H and high_M and low_R:
        return "malignant_reservoir"

    if low_H and low_M and high_R and high_S:
        return "stable_recoverable"

    if low_H and low_M and high_R:
        return "recoverable_reference"

    if (abs(row["H"] - med["H"]) < 0.3) and (abs(row["M"] - med["M"]) < 0.3):
        return "transitional_intermediate"

    return "mixed_sector"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    tables_dir = Path(cfg["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    state_path = processed_dir / "state_table_refined.csv"
    if not state_path.exists():
        raise FileNotFoundError(f"Missing refined state table: {state_path}")

    df = read_table(state_path).copy()

    med = {
        "H": float(df["H"].median()),
        "S": float(df["S"].median()),
        "M": float(df["M"].median()),
        "R": float(df["R"].median()),
    }

    df["branch_sector"] = df.apply(assign_sector, axis=1, med=med)

    out_path = tables_dir / "branch_sector_assignments.csv"
    write_table(df, out_path)

    sector_summary = (
        df.groupby(["population_label", "branch_sector"])
        .size()
        .reset_index(name="n_cells")
        .sort_values(["population_label", "n_cells"], ascending=[True, False])
    )
    sector_summary_path = tables_dir / "branch_sector_summary.csv"
    write_table(sector_summary, sector_summary_path)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "medians_used": med,
        "outputs": {
            "assignments": str(out_path),
            "summary": str(sector_summary_path),
        },
        "sector_counts": df["branch_sector"].value_counts().to_dict(),
    }

    with open(tables_dir / "branch_sector_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
