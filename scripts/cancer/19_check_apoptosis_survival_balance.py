import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

from src.common.score_utils import (
    expression_wide_to_cell_by_gene,
    attach_gene_symbols_to_cell_by_gene,
    score_gene_set_mean_z,
    gene_set_match_report,
)
from src.common.io import read_table, write_table


def load_yaml(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--gene-set-config", default="configs/overlays/gene_sets_human_basic.yaml")
    args = parser.parse_args()

    cfg = load_yaml(args.config)
    gs = load_yaml(args.gene_set_config)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    tables_dir = Path(cfg["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    expr = read_table(processed_dir / "expression_matrix_wide.parquet")
    state = read_table(processed_dir / "state_table.csv")
    gene_map = read_table(processed_dir / "gene_map.csv")

    x = expression_wide_to_cell_by_gene(expr, gene_col="ensembl_id")
    x = attach_gene_symbols_to_cell_by_gene(
        x,
        gene_map,
        ensembl_col="ensembl_id",
        symbol_col="gene_symbol",
    )

    surv_genes = gs["apoptosis_survival"]["survival"]
    apop_genes = gs["apoptosis_survival"]["apoptosis"]

    survival = score_gene_set_mean_z(x, surv_genes)
    apoptosis = score_gene_set_mean_z(x, apop_genes)
    balance = survival - apoptosis

    overlay = pd.DataFrame(
        {
            "cell_id": x.index,
            "survival_score": survival.values,
            "apoptosis_score": apoptosis.values,
            "survival_minus_apoptosis_score": balance.values,
        }
    )

    merged = state.merge(overlay, on="cell_id", how="left")

    overlay_path = tables_dir / "apoptosis_survival_overlay_table.csv"
    write_table(merged, overlay_path)

    group_summary = (
        merged.groupby("population_label")[
            ["survival_score", "apoptosis_score", "survival_minus_apoptosis_score", "H", "S", "M", "R"]
        ]
        .mean()
        .round(4)
        .reset_index()
    )
    group_summary_path = tables_dir / "apoptosis_survival_overlay_group_summary.csv"
    write_table(group_summary, group_summary_path)

    corr = (
        merged[
            ["survival_score", "apoptosis_score", "survival_minus_apoptosis_score", "H", "S", "M", "R"]
        ]
        .corr()
        .round(4)
        .reset_index()
    )
    corr_path = tables_dir / "apoptosis_survival_overlay_correlations.csv"
    write_table(corr, corr_path)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "survival_match_report": gene_set_match_report(x, surv_genes),
        "apoptosis_match_report": gene_set_match_report(x, apop_genes),
        "outputs": {
            "overlay_table": str(overlay_path),
            "group_summary": str(group_summary_path),
            "correlations": str(corr_path),
        },
    }

    with open(tables_dir / "apoptosis_survival_overlay_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
