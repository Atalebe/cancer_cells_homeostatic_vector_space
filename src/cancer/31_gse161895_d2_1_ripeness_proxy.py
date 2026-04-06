#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml
from scipy.stats import mannwhitneyu


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def rank01(s: pd.Series) -> pd.Series:
    return s.rank(method="average", pct=True)


def main() -> None:
    cfg = read_config()
    proc_dir = REPO_ROOT / cfg["processed_dir"]
    d2_dir = REPO_ROOT / cfg["results_dir"] / "d2_only_reanalysis"
    out_dir = REPO_ROOT / cfg["results_dir"] / "d2_1_ripeness_proxy"
    out_dir.mkdir(parents=True, exist_ok=True)

    state = pd.read_parquet(proc_dir / "state_table.parquet")[["cell_id", "H", "S", "M", "R", "phi"]]
    meta = pd.read_parquet(proc_dir / "cell_metadata_registry_aligned.parquet")[["cell_id", "treatment_state"]]
    d2 = pd.read_parquet(d2_dir / "d2_state_domains.parquet")[["cell_id", "d2_subdomain"]]

    df = state.merge(meta, on="cell_id", how="left").merge(d2, on="cell_id", how="left")
    df = df.loc[
        (df["d2_subdomain"] == "D2_1") &
        (df["treatment_state"].isin(["treated", "untreated"]))
    ].copy()

    df["H_rank"] = rank01(df["H"])
    df["M_rank"] = rank01(df["M"])
    df["phi_rank"] = rank01(df["phi"])
    df["inv_R_rank"] = 1.0 - rank01(df["R"])

    df["ripeness_proxy"] = (
        df["H_rank"] + df["M_rank"] + df["phi_rank"] + df["inv_R_rank"]
    ) / 4.0

    df.to_csv(out_dir / "d2_1_ripeness_proxy_cells.csv", index=False)

    a = df.loc[df["treatment_state"] == "treated", "ripeness_proxy"]
    b = df.loc[df["treatment_state"] == "untreated", "ripeness_proxy"]
    u, p = mannwhitneyu(a, b, alternative="two-sided")

    summary = {
        "definition": "rank-based internal proxy = mean(rank(H), rank(M), rank(phi), 1-rank(R)) within D2_1",
        "treated_n": int(len(a)),
        "untreated_n": int(len(b)),
        "median_treated": float(a.median()),
        "median_untreated": float(b.median()),
        "delta_median_treated_minus_untreated": float(a.median() - b.median()),
        "u_statistic": float(u),
        "p_value": float(p),
        "note": "This is an internal state-ripeness sensitivity analysis, not chronological age."
    }

    with open(out_dir / "d2_1_ripeness_proxy_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
