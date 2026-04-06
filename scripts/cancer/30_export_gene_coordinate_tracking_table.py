import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

from src.common.io import read_table, write_table


def strongest_coordinate(row) -> str:
    vals = {
        "H": abs(row.get("corr_H", 0.0)),
        "S": abs(row.get("corr_S", 0.0)),
        "M": abs(row.get("corr_M", 0.0)),
        "R": abs(row.get("corr_R", 0.0)),
    }
    return max(vals, key=vals.get)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--min-lists", type=int, default=3)
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    tables_dir = Path(cfg["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    core = read_table(tables_dir / "irreversible_core_candidate_top_ranked.csv")
    corr = read_table(tables_dir / "gene_coordinate_correlations.csv")

    core_sub = core[core["n_lists_present"] >= args.min_lists].head(args.top_n).copy()
    corr = corr.rename(columns={"gene": "gene_symbol"})
    core_sub = core_sub.rename(columns={"gene": "gene_symbol"})

    merged = core_sub.merge(corr, on="gene_symbol", how="left")
    merged["strongest_coordinate"] = merged.apply(strongest_coordinate, axis=1)

    out_cols = [
        "gene_symbol",
        "n_lists_present",
        "present_in_lists",
        "mean_rank_across_present_lists",
        "corr_H",
        "corr_S",
        "corr_M",
        "corr_R",
        "strongest_coordinate",
    ]
    out_df = merged[out_cols].copy()

    out_path = tables_dir / "candidate_core_gene_coordinate_tracking.csv"
    write_table(out_df, out_path)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "min_lists": int(args.min_lists),
        "top_n": int(args.top_n),
        "n_rows": int(out_df.shape[0]),
        "output_file": str(out_path),
    }

    with open(tables_dir / "candidate_core_gene_coordinate_tracking_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
