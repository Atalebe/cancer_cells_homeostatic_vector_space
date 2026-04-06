import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from src.common.io import read_table, write_table
from src.common.score_utils import (
    expression_wide_to_cell_by_gene,
    attach_gene_symbols_to_cell_by_gene,
)


def robust_center_rows(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for idx in out.index:
        row = out.loc[idx].astype(float)
        med = np.median(row.values)
        mad = np.median(np.abs(row.values - med))
        if mad == 0 or np.isnan(mad):
            out.loc[idx] = row.values - med
        else:
            out.loc[idx] = (row.values - med) / (1.4826 * mad)
    return out


def plot_heatmap(df: pd.DataFrame, out_path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(10, max(6, 0.35 * len(df.index))))
    im = ax.imshow(df.values, aspect="auto")
    ax.set_xticks(range(len(df.columns)))
    ax.set_xticklabels(df.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(df.index)))
    ax.set_yticklabels(df.index)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--min-lists", type=int, default=3)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    tables_dir = Path(cfg["paths"]["tables_dir"])
    figures_dir = Path(cfg["paths"]["figures_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    expr_path = processed_dir / "expression_matrix_wide.parquet"
    gene_map_path = processed_dir / "gene_map.csv"
    refined_state_path = processed_dir / "state_table_refined.csv"
    branch_path = tables_dir / "branch_sector_assignments.csv"
    core_path = tables_dir / "irreversible_core_candidate_top_ranked.csv"

    if not expr_path.exists():
        raise FileNotFoundError(f"Missing expression matrix: {expr_path}")
    if not gene_map_path.exists():
        raise FileNotFoundError(f"Missing gene map: {gene_map_path}")
    if not refined_state_path.exists():
        raise FileNotFoundError(f"Missing refined state table: {refined_state_path}")
    if not branch_path.exists():
        raise FileNotFoundError(f"Missing branch sector assignments: {branch_path}")
    if not core_path.exists():
        raise FileNotFoundError(f"Missing candidate core table: {core_path}")

    expr = read_table(expr_path)
    gene_map = read_table(gene_map_path)
    state = read_table(refined_state_path)
    branch = read_table(branch_path)
    core = read_table(core_path)

    x = expression_wide_to_cell_by_gene(expr, gene_col="ensembl_id")
    x = attach_gene_symbols_to_cell_by_gene(
        x,
        gene_map,
        ensembl_col="ensembl_id",
        symbol_col="gene_symbol",
    )
    x_log = np.log1p(x)

    meta = state.merge(branch[["cell_id", "branch_sector"]], on="cell_id", how="left").set_index("cell_id")
    x_log = x_log.loc[meta.index]

    core_sub = core[core["n_lists_present"] >= args.min_lists].copy()
    core_sub = core_sub.head(args.top_n).copy()
    selected_genes = [g for g in core_sub["gene"].tolist() if g in x_log.columns]

    if len(selected_genes) == 0:
        raise ValueError("No selected candidate genes matched expression matrix columns.")

    # population means
    pop_mat = (
        x_log[selected_genes]
        .groupby(meta["population_label"])
        .mean()
        .T
    )
    pop_mat_z = robust_center_rows(pop_mat)

    pop_csv = tables_dir / "candidate_core_heatmap_population_means.csv"
    pop_z_csv = tables_dir / "candidate_core_heatmap_population_means_rowz.csv"
    write_table(pop_mat.reset_index().rename(columns={"index": "gene"}), pop_csv)
    write_table(pop_mat_z.reset_index().rename(columns={"index": "gene"}), pop_z_csv)

    # sector means
    sector_mat = (
        x_log[selected_genes]
        .groupby(meta["branch_sector"])
        .mean()
        .T
    )
    sector_mat_z = robust_center_rows(sector_mat)

    sector_csv = tables_dir / "candidate_core_heatmap_sector_means.csv"
    sector_z_csv = tables_dir / "candidate_core_heatmap_sector_means_rowz.csv"
    write_table(sector_mat.reset_index().rename(columns={"index": "gene"}), sector_csv)
    write_table(sector_mat_z.reset_index().rename(columns={"index": "gene"}), sector_z_csv)

    # plots
    pop_fig = figures_dir / "candidate_core_heatmap_by_population.png"
    pop_z_fig = figures_dir / "candidate_core_heatmap_by_population_rowz.png"
    sector_fig = figures_dir / "candidate_core_heatmap_by_sector.png"
    sector_z_fig = figures_dir / "candidate_core_heatmap_by_sector_rowz.png"

    plot_heatmap(pop_mat, pop_fig, "Candidate core genes by population")
    plot_heatmap(pop_mat_z, pop_z_fig, "Candidate core genes by population, row-z")
    plot_heatmap(sector_mat, sector_fig, "Candidate core genes by branch sector")
    plot_heatmap(sector_mat_z, sector_z_fig, "Candidate core genes by branch sector, row-z")

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "top_n": int(args.top_n),
        "min_lists": int(args.min_lists),
        "n_selected_genes": int(len(selected_genes)),
        "selected_genes": selected_genes,
        "outputs": {
            "population_means_csv": str(pop_csv),
            "population_means_rowz_csv": str(pop_z_csv),
            "sector_means_csv": str(sector_csv),
            "sector_means_rowz_csv": str(sector_z_csv),
            "population_heatmap": str(pop_fig),
            "population_heatmap_rowz": str(pop_z_fig),
            "sector_heatmap": str(sector_fig),
            "sector_heatmap_rowz": str(sector_z_fig),
        },
    }

    with open(tables_dir / "candidate_core_heatmap_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
