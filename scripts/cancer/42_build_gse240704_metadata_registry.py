import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

from src.common.io import read_table, write_table


def infer_population_label(cell_id: str) -> str:
    low = cell_id.lower()
    if "stem" in low or "malignant" in low or "tumor" in low:
        return "putative_malignant"
    if "normal" in low or "ctrl" in low or "control" in low:
        return "putative_reference"
    return "unknown"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--matrix-kind", choices=["normalized", "non_normalized"], default="normalized")
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    metadata_dir = Path(cfg["paths"]["metadata_dir"])
    metadata_dir.mkdir(parents=True, exist_ok=True)

    expr = read_table(processed_dir / f"gse240704_expression_{args.matrix_kind}.parquet")
    cell_cols = [c for c in expr.columns if c != "gene_id"]

    meta = pd.DataFrame({"cell_id": cell_cols})
    meta["population_label"] = meta["cell_id"].astype(str).apply(infer_population_label)
    meta["dataset_id"] = cfg["dataset"]["dataset_id"]
    meta["accession"] = cfg["dataset"]["accession"]

    out_csv = metadata_dir / f"gse240704_metadata_registry_{args.matrix_kind}.csv"
    write_table(meta, out_csv)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "matrix_kind": args.matrix_kind,
        "n_cells": int(meta.shape[0]),
        "population_counts": meta["population_label"].value_counts().to_dict(),
        "output_file": str(out_csv),
    }

    with open(metadata_dir / f"gse240704_metadata_registry_{args.matrix_kind}_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
