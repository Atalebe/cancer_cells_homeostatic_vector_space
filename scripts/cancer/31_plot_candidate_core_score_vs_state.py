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
    fig, ax = plt.subplots(figsize=(8, 6))
    grouped = list(df.groupby(group_col))
    groups = [sub[value_col].values for _, sub in grouped]
    labels = [label for label, _ in grouped]
    ax.boxplot(groups, tick_labels=labels)
    ax.set_ylabel(value_col)
    ax.set_title(title)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    tables_dir = Path(cfg["paths"]["tables_dir"])
    figures_dir = Path(cfg["paths"]["figures_dir"])
    figures_dir.mkdir(parents=True, exist_ok=True)

    score_path = tables_dir / "candidate_core_score_table.csv"
    if not score_path.exists():
        raise FileNotFoundError(f"Missing candidate core score table: {score_path}")

    df = read_table(score_path)

    # population-colored scatter plots
    for coord in ["H", "S", "M", "R"]:
        scatter_plot(
            df,
            x="candidate_core_score",
            y=coord,
            color_col="population_label",
            out_path=figures_dir / f"candidate_core_score_vs_{coord}_by_population.png",
            title=f"Candidate core score vs {coord}, by population",
        )

    # sector-colored scatter plots
    for coord in ["H", "S", "M", "R"]:
        scatter_plot(
            df,
            x="candidate_core_score",
            y=coord,
            color_col="branch_sector",
            out_path=figures_dir / f"candidate_core_score_vs_{coord}_by_sector.png",
            title=f"Candidate core score vs {coord}, by sector",
        )

    # boxplots
    box_plot(
        df,
        value_col="candidate_core_score",
        group_col="population_label",
        out_path=figures_dir / "candidate_core_score_boxplot_by_population.png",
        title="Candidate core score by population",
    )

    box_plot(
        df,
        value_col="candidate_core_score",
        group_col="branch_sector",
        out_path=figures_dir / "candidate_core_score_boxplot_by_sector.png",
        title="Candidate core score by branch sector",
    )

    print(f"[ok] wrote candidate-core score figures to {figures_dir}")


if __name__ == "__main__":
    main()
