#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import yaml


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def resolve_state_domain_column(df: pd.DataFrame) -> str:
    candidates = ["state_domain", "domain", "cluster", "label"]
    for col in candidates:
        if col in df.columns:
            return col
    raise ValueError(f"Could not find state-domain column in columns: {list(df.columns)}")


def main() -> None:
    cfg = read_config()
    proc_dir = REPO_ROOT / cfg["processed_dir"]
    results_dir = REPO_ROOT / cfg["results_dir"]
    d2_dir = results_dir / "d2_only_reanalysis"
    state_domains_path = results_dir / "state_domains" / "state_domains.parquet"

    out_dir = results_dir / "hrsm_profile_visual_compare"
    plot_dir = out_dir / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    state = pd.read_parquet(proc_dir / "state_table.parquet")[["cell_id", "H", "S", "M", "R", "phi"]]

    state_domains = pd.read_parquet(state_domains_path)
    state_domain_col = resolve_state_domain_column(state_domains)
    state_domains = state_domains.rename(columns={state_domain_col: "state_domain"})
    state_domains = state_domains[["cell_id", "state_domain"]]

    meta = pd.read_parquet(proc_dir / "cell_metadata_registry_aligned.parquet")[
        ["cell_id", "patient_or_donor_id", "source_class", "treatment_state"]
    ]
    d2 = pd.read_parquet(d2_dir / "d2_state_domains.parquet")[["cell_id", "d2_subdomain"]]

    df = (
        state
        .merge(state_domains, on="cell_id", how="left")
        .merge(meta, on="cell_id", how="left")
        .merge(d2, on="cell_id", how="left")
    )

    groups = [
        ("D2_1 treated", (df["d2_subdomain"] == "D2_1") & (df["treatment_state"] == "treated")),
        ("D2_1 untreated", (df["d2_subdomain"] == "D2_1") & (df["treatment_state"] == "untreated")),
        ("D2_2 immune-like", df["d2_subdomain"] == "D2_2"),
        ("Normal donors", df["source_class"] == "normal_donor"),
    ]

    frames = []
    summary_rows = []

    for label, mask in groups:
        sub = df.loc[mask].copy()
        if sub.empty:
            continue
        sub["comparison_group"] = label
        frames.append(sub)

        row = {
            "comparison_group": label,
            "n_cells": int(len(sub)),
        }
        for metric in ["H", "S", "M", "R", "phi"]:
            row[f"{metric}_median"] = float(sub[metric].median())
            row[f"{metric}_mean"] = float(sub[metric].mean())
        summary_rows.append(row)

    if not frames:
        raise ValueError("No comparison groups contained cells.")

    comp = pd.concat(frames, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_dir / "hrsm_profile_group_summary.csv", index=False)
    comp.to_csv(out_dir / "hrsm_profile_cells_input.csv", index=False)

    metrics = ["H", "S", "M", "R", "phi"]
    long_df = comp.melt(
        id_vars=["cell_id", "comparison_group"],
        value_vars=metrics,
        var_name="metric",
        value_name="value",
    )
    long_df.to_csv(out_dir / "hrsm_profile_long.csv", index=False)

    # Grouped boxplot across all metrics
    fig, ax = plt.subplots(figsize=(14, 7))

    group_names = summary["comparison_group"].tolist()
    n_groups = len(group_names)
    positions = []
    data_series = []
    tick_positions = []
    tick_labels = []

    width = 0.18
    offsets = []
    if n_groups == 1:
        offsets = [0.0]
    else:
        start = -width * (n_groups - 1) / 2.0
        offsets = [start + i * width for i in range(n_groups)]

    for i, metric in enumerate(metrics):
        tick_positions.append(i + 1)
        tick_labels.append(metric)
        for j, group_name in enumerate(group_names):
            vals = long_df.loc[
                (long_df["metric"] == metric) &
                (long_df["comparison_group"] == group_name),
                "value"
            ].dropna().to_numpy()
            positions.append((i + 1) + offsets[j])
            data_series.append(vals)

    ax.boxplot(
        data_series,
        positions=positions,
        widths=width * 0.9,
        patch_artist=False,
        showfliers=False,
        manage_ticks=False,
    )

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.set_ylabel("HRSM value")
    ax.set_title("HRSM profiles across treated, untreated, immune-like, and donor groups")

    legend_labels = []
    for j, group_name in enumerate(group_names):
        legend_labels.append(f"{group_name}")
    ax.legend(legend_labels, loc="best")

    plt.tight_layout()
    plt.savefig(plot_dir / "hrsm_profile_boxplot_by_group.png", dpi=200)
    plt.close()

    # Per-metric boxplots
    for metric in metrics:
        fig, ax = plt.subplots(figsize=(9, 5))
        data = []
        labels = []
        for group_name in group_names:
            vals = comp.loc[comp["comparison_group"] == group_name, metric].dropna().to_numpy()
            data.append(vals)
            labels.append(group_name)

        ax.boxplot(data, labels=labels, showfliers=False)
        ax.set_ylabel(metric)
        ax.set_title(f"{metric} by comparison group")
        plt.xticks(rotation=20, ha="right")
        plt.tight_layout()
        plt.savefig(plot_dir / f"{metric.lower()}_by_comparison_group.png", dpi=200)
        plt.close()

    with open(out_dir / "hrsm_profile_visual_compare_summary.json", "w", encoding="utf-8") as f:
        json.dump(
            {
                "state_domain_source": str(state_domains_path),
                "groups_written": group_names,
                "n_groups": int(len(group_names)),
                "outputs": {
                    "summary_csv": str(out_dir / "hrsm_profile_group_summary.csv"),
                    "long_csv": str(out_dir / "hrsm_profile_long.csv"),
                    "boxplot": str(plot_dir / "hrsm_profile_boxplot_by_group.png"),
                },
            },
            f,
            indent=2,
        )

    print(summary)


if __name__ == "__main__":
    main()
