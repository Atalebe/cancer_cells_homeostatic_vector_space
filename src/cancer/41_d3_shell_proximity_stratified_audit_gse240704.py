#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu


def median_or_nan(s: pd.Series) -> float:
    s = pd.to_numeric(s, errors="coerce")
    if s.dropna().empty:
        return np.nan
    return float(s.median())


def mwu(a: pd.Series, b: pd.Series):
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) == 0 or len(b) == 0:
        return np.nan, np.nan
    stat, p = mannwhitneyu(a, b, alternative="two-sided")
    return float(stat), float(p)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-level-csv", required=True)
    ap.add_argument("--proximity-csv", required=True)
    ap.add_argument("--outcsv", required=True)
    ap.add_argument("--outjson", required=True)
    args = ap.parse_args()

    sample_df = pd.read_csv(args.sample_level_csv)
    prox_df = pd.read_csv(args.proximity_csv)

    if "sample_id" not in sample_df.columns:
        raise ValueError("sample-level csv must contain sample_id")
    if "sample_id" not in prox_df.columns:
        raise ValueError("proximity csv must contain sample_id")

    merged = sample_df.merge(prox_df, on="sample_id", how="left", suffixes=("", "_prox"))

    required = ["state_domain", "phi", "H", "S", "M", "R"]
    missing = [c for c in required if c not in merged.columns]
    if missing:
        raise ValueError(f"Missing required columns after merge: {missing}")

    distance_col = None
    distance_candidates = [
        "min_hrsm_distance_to_shell",
        "min_distance_to_stable_shell",
        "min_distance_to_shell",
        "distance_to_stable_shell",
    ]
    for c in distance_candidates:
        if c in merged.columns:
            distance_col = c
            break
    if distance_col is None:
        raise ValueError(
            "Could not find shell distance column. Checked: "
            + ", ".join(distance_candidates)
        )

    d3 = merged[merged["state_domain"] == "D3"].copy()
    if d3.empty:
        raise ValueError("No D3 rows found")

    d3[distance_col] = pd.to_numeric(d3[distance_col], errors="coerce")
    d3 = d3[d3[distance_col].notna()].copy()
    if d3.empty:
        raise ValueError("D3 rows found, but all shell distances are missing")

    cutoff = float(d3[distance_col].median())
    d3["d3_shell_proximity_group"] = np.where(
        d3[distance_col] <= cutoff,
        "D3_shell_near",
        "D3_shell_far",
    )

    cond_col = "placeholder_condition_cur" if "placeholder_condition_cur" in d3.columns else None
    if cond_col is not None:
        d3["_cond"] = d3[cond_col].fillna("nan").astype(str)
    else:
        d3["_cond"] = "nan"

    near = d3[d3["d3_shell_proximity_group"] == "D3_shell_near"].copy()
    far = d3[d3["d3_shell_proximity_group"] == "D3_shell_far"].copy()

    rows = []

    for axis in ["H", "S", "M", "R", "phi"]:
        stat, p = mwu(near[axis], far[axis])
        rows.append({
            "comparison": f"D3_shell_near_vs_far_{axis}",
            "axis": axis,
            "n_near": int(len(near)),
            "n_far": int(len(far)),
            "median_near": median_or_nan(near[axis]),
            "median_far": median_or_nan(far[axis]),
            "delta_median_near_minus_far": median_or_nan(near[axis]) - median_or_nan(far[axis]),
            "u_statistic": stat,
            "p_value": p,
        })

    for grp_name, sub in d3.groupby("d3_shell_proximity_group", dropna=False):
        rows.append({
            "comparison": f"{grp_name}_condition_counts",
            "axis": "",
            "n_near": int(len(sub)) if grp_name == "D3_shell_near" else np.nan,
            "n_far": int(len(sub)) if grp_name == "D3_shell_far" else np.nan,
            "median_near": np.nan,
            "median_far": np.nan,
            "delta_median_near_minus_far": np.nan,
            "u_statistic": np.nan,
            "p_value": np.nan,
            "condition_counts_json": json.dumps(sub["_cond"].value_counts().to_dict()),
        })

    meth_near = int((near["_cond"] == "mgmt_methylated").sum())
    nonmeth_near = int((near["_cond"] == "mgmt_non_methylated").sum())
    meth_far = int((far["_cond"] == "mgmt_methylated").sum())
    nonmeth_far = int((far["_cond"] == "mgmt_non_methylated").sum())

    if (meth_near + nonmeth_near) > 0 and (meth_far + nonmeth_far) > 0:
        table = [[meth_near, nonmeth_near], [meth_far, nonmeth_far]]
        odds, p = fisher_exact(table)
        rows.append({
            "comparison": "D3_shell_near_vs_far_mgmt_enrichment",
            "axis": "condition",
            "n_near": int(len(near)),
            "n_far": int(len(far)),
            "median_near": np.nan,
            "median_far": np.nan,
            "delta_median_near_minus_far": np.nan,
            "u_statistic": np.nan,
            "p_value": float(p),
            "fisher_odds_ratio": float(odds),
            "meth_near": meth_near,
            "nonmeth_near": nonmeth_near,
            "meth_far": meth_far,
            "nonmeth_far": nonmeth_far,
        })

    out = pd.DataFrame(rows)
    Path(args.outcsv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.outcsv, index=False)

    summary = {
        "n_total_rows_sample_level": int(len(sample_df)),
        "n_total_rows_proximity": int(len(prox_df)),
        "n_d3_after_merge": int(len(d3)),
        "distance_column_used": distance_col,
        "d3_shell_distance_cutoff_median": cutoff,
        "n_d3_shell_near": int(len(near)),
        "n_d3_shell_far": int(len(far)),
        "condition_column_used": cond_col,
    }
    with open(args.outjson, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[ok] wrote D3 shell proximity audit: {args.outcsv}")
    print(f"[ok] wrote summary: {args.outjson}")
    print("[info] summary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
