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


def mwu_test(a: pd.Series, b: pd.Series) -> dict:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) == 0 or len(b) == 0:
        return {
            "n_a": len(a),
            "n_b": len(b),
            "median_a": None,
            "median_b": None,
            "delta_median_a_minus_b": None,
            "u_statistic": None,
            "p_value": None,
        }
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(a.median()),
        "median_b": float(b.median()),
        "delta_median_a_minus_b": float(a.median() - b.median()),
        "u_statistic": float(u),
        "p_value": float(p),
    }


def main() -> None:
    cfg = read_config()
    proc_dir = REPO_ROOT / cfg["processed_dir"]
    state_dir = REPO_ROOT / cfg["results_dir"] / "state_domains"
    out_dir = REPO_ROOT / cfg["results_dir"] / "domain_hrsm_tests"
    out_dir.mkdir(parents=True, exist_ok=True)

    state = pd.read_parquet(proc_dir / "state_table.parquet")
    domains = pd.read_parquet(state_dir / "state_domains.parquet")[["cell_id", "state_domain"]]
    df = state.merge(domains, on="cell_id", how="inner")

    rows = []
    for metric in ["H", "S", "M", "R", "phi"]:
        a = df.loc[df["state_domain"] == "D1", metric]
        b = df.loc[df["state_domain"] == "D2", metric]
        res = mwu_test(a, b)
        res["metric"] = metric
        res["group_a"] = "D1"
        res["group_b"] = "D2"
        rows.append(res)

    out = pd.DataFrame(rows)
    out = out[[
        "metric", "group_a", "group_b",
        "n_a", "n_b",
        "median_a", "median_b", "delta_median_a_minus_b",
        "u_statistic", "p_value"
    ]]
    out.to_csv(out_dir / "d1_vs_d2_hrsm_tests.csv", index=False)

    with open(out_dir / "d1_vs_d2_hrsm_tests.json", "w", encoding="utf-8") as f:
        json.dump(out.to_dict(orient="records"), f, indent=2)

    print(out)


if __name__ == "__main__":
    main()
