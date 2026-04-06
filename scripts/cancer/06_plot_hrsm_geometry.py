import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import yaml

from src.common.io import read_table


def scatter_plot(df, x, y, color_col, out_path, title):
    fig, ax = plt.subplots(figsize=(7, 6))
    for label, sub in df.groupby(color_col):
        ax.scatter(sub[x], sub[y], label=label, alpha=0.8)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def box_plot(df, value_col, group_col, out_path, title):
    fig, ax = plt.subplots(figsize=(7, 6))
    groups = [sub[value_col].values for _, sub in df.groupby(group_col)]
    labels = [label for label, _ in df.groupby(group_col)]
    ax.boxplot(groups, tick_labels=labels)
    ax.set_ylabel(value_col)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    figures_dir = Path(cfg["paths"]["figures_dir"])
    figures_dir.mkdir(parents=True, exist_ok=True)

    state_path = processed_dir / "state_table_refined.csv"
    if not state_path.exists():
        raise FileNotFoundError(f"Missing refined state table: {state_path}")

    df = read_table(state_path)

    scatter_plot(df, "H", "R", "population_label", figures_dir / "hrsm_H_vs_R_refined.png", "Refined HRSM: H vs R")
    scatter_plot(df, "M", "R", "population_label", figures_dir / "hrsm_M_vs_R_refined.png", "Refined HRSM: M vs R")
    scatter_plot(df, "S", "R", "population_label", figures_dir / "hrsm_S_vs_R_refined.png", "Refined HRSM: S vs R")
    scatter_plot(df, "pc1", "pc2", "population_label", figures_dir / "pca_pc1_vs_pc2_by_population.png", "PCA: PC1 vs PC2")

    for coord in ["H", "S", "M", "R"]:
        box_plot(
            df,
            value_col=coord,
            group_col="population_label",
            out_path=figures_dir / f"boxplot_{coord}_by_population.png",
            title=f"{coord} by population"
        )

    print(f"[ok] wrote figures to {figures_dir}")


if __name__ == "__main__":
    main()
