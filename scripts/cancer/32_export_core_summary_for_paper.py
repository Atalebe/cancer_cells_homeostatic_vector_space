import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

from src.common.io import read_table, write_table


def assign_theme(gene: str) -> str:
    themes = {
        "FN1": "matrix_adhesion",
        "COL6A1": "matrix_adhesion",
        "LGALS3BP": "matrix_adhesion",
        "MARCKSL1": "cytoskeletal_plasticity",
        "GSN": "cytoskeletal_plasticity",
        "NET1": "cytoskeletal_plasticity",
        "ARL4C": "membrane_plasticity",
        "ABCC3": "transport_adaptation",
        "SLCO4A1": "transport_adaptation",
        "GALNT10": "glycosylation_remodeling",
        "GALNT5": "glycosylation_remodeling",
        "DDX17": "rna_processing_regulation",
        "MAGI1": "signaling_scaffold",
        "ABR": "signaling_regulation",
        "TMEM131": "membrane_trafficking",
        "ESYT2": "membrane_contact_lipid_transfer",
    }
    return themes.get(gene, "other")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--min-lists", type=int, default=3)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    tables_dir = Path(cfg["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    interp_path = tables_dir / "candidate_core_interpretation_table.csv"
    track_path = tables_dir / "candidate_core_gene_coordinate_tracking.csv"

    if not interp_path.exists():
        raise FileNotFoundError(f"Missing interpretation table: {interp_path}")
    if not track_path.exists():
        raise FileNotFoundError(f"Missing tracking table: {track_path}")

    interp = read_table(interp_path)
    track = read_table(track_path)

    df = interp.merge(
        track[["gene_symbol", "strongest_coordinate"]],
        on="gene_symbol",
        how="left",
    )

    df = df[df["n_lists_present"] >= args.min_lists].copy().head(args.top_n)
    df["biological_theme"] = df["gene_symbol"].astype(str).apply(assign_theme)

    out_cols = [
        "gene_symbol",
        "n_lists_present",
        "present_in_lists",
        "strongest_coordinate",
        "corr_H",
        "corr_M",
        "corr_R",
        "malignant_reservoir_logfc",
        "unstable_committed_logfc",
        "stable_recoverable_logfc",
        "PKH26_High_vs_G1_logfc",
        "biological_theme",
        "interpretation_tags",
    ]
    out_df = df[out_cols].copy()

    out_path = tables_dir / "candidate_core_summary_for_paper.csv"
    write_table(out_df, out_path)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "top_n": int(args.top_n),
        "min_lists": int(args.min_lists),
        "n_rows": int(out_df.shape[0]),
        "output_file": str(out_path),
    }

    with open(tables_dir / "candidate_core_summary_for_paper_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
