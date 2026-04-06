import argparse
import json
from pathlib import Path

import pandas as pd
import yaml


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    raw_dir = Path(cfg["paths"]["raw_dir"])
    metadata_dir = Path(cfg["paths"]["metadata_dir"])
    metadata_dir.mkdir(parents=True, exist_ok=True)

    xlsx_path = raw_dir / "GSE124989_Reads.xlsx"
    if not xlsx_path.exists():
        raise FileNotFoundError(f"Missing workbook: {xlsx_path}")

    df0 = pd.read_excel(xlsx_path, sheet_name="Blad1", nrows=5)
    columns = list(df0.columns)

    if len(columns) < 2:
        raise ValueError("Workbook does not look like a gene x cell matrix.")

    cell_ids = [str(c) for c in columns[1:]]

    rows = []
    for cell_id in cell_ids:
        prefix = cell_id.split("-")[0]

        if prefix == "G1":
            population_label = "G1"
            group_order = 0
        elif "HIGH" in cell_id.upper():
            population_label = "PKH26_High"
            group_order = 2
        elif "LOW" in cell_id.upper():
            population_label = "PKH26_Low"
            group_order = 1
        else:
            population_label = "UNKNOWN"
            group_order = -1

        rows.append(
            {
                "cell_id": cell_id,
                "raw_prefix": prefix,
                "population_label": population_label,
                "group_order": group_order,
                "dataset_id": cfg["dataset"]["dataset_id"],
                "accession": cfg["dataset"]["accession"],
            }
        )

    meta = pd.DataFrame(rows)

    meta_path = metadata_dir / "cell_metadata_registry.csv"
    meta.to_csv(meta_path, index=False)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "n_cells": int(meta.shape[0]),
        "population_counts": meta["population_label"].value_counts(dropna=False).to_dict(),
        "unknown_labels": int((meta["population_label"] == "UNKNOWN").sum()),
        "output_file": str(meta_path),
    }

    with open(metadata_dir / "cell_metadata_registry_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
