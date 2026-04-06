import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

from src.common.io import read_table, write_table


def assign_sector(row, h_cut, s_cut, m_cut, r_cut):
    high_H = row["H"] >= h_cut
    high_M = row["M"] >= m_cut
    high_S = row["S"] >= s_cut
    high_R = row["R"] >= r_cut

    low_H = row["H"] < -h_cut
    low_M = row["M"] < -m_cut
    low_S = row["S"] < -s_cut
    low_R = row["R"] < -r_cut

    near_H = abs(row["H"]) < h_cut
    near_M = abs(row["M"]) < m_cut
    near_R = abs(row["R"]) < r_cut

    if high_H and high_M and low_R and low_S:
        return "unstable_committed_malignant"
    if high_H and high_M and low_R:
        return "malignant_reservoir"
    if low_H and low_M and high_R and high_S:
        return "stable_recoverable"
    if low_H and low_M and high_R:
        return "recoverable_reference"
    if near_H and near_M and near_R:
        return "transitional_intermediate"
    return "mixed_sector"


def run_scheme(df: pd.DataFrame, scheme_name: str, cut: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = df.copy()
    out["threshold_scheme"] = scheme_name
    out["branch_sector"] = out.apply(
        assign_sector,
        axis=1,
        h_cut=cut,
        s_cut=cut,
        m_cut=cut,
        r_cut=cut,
    )

    summary = (
        out.groupby(["population_label", "branch_sector"])
        .size()
        .reset_index(name="n_cells")
        .sort_values(["population_label", "n_cells"], ascending=[True, False])
    )
    summary["threshold_scheme"] = scheme_name
    return out, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    tables_dir = Path(cfg["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    state = read_table(processed_dir / "state_table_refined.csv")

    schemes = {
        "median_like_0.0": 0.0,
        "moderate_0.5": 0.5,
        "strict_0.75": 0.75,
    }

    all_assignments = []
    all_summaries = []

    for name, cut in schemes.items():
        assigned, summary = run_scheme(state, name, cut)
        all_assignments.append(assigned)
        all_summaries.append(summary)

    assign_df = pd.concat(all_assignments, ignore_index=True)
    summary_df = pd.concat(all_summaries, ignore_index=True)

    assign_path = tables_dir / "branch_threshold_robustness_assignments.csv"
    summary_path = tables_dir / "branch_threshold_robustness_summary.csv"
    write_table(assign_df, assign_path)
    write_table(summary_df, summary_path)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "schemes": schemes,
        "output_assignments": str(assign_path),
        "output_summary": str(summary_path),
    }

    with open(tables_dir / "branch_threshold_robustness_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
