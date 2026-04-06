#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Doctor-style audit of D3 as an outlier branch.")
    p.add_argument("--selected-probe-matrix", required=True)
    p.add_argument("--merged-curated-csv", required=True)
    p.add_argument("--stable-shell-csv", required=True)
    p.add_argument("--neighbor-csv", required=True, help="stable_shell_k_nearest_neighbors.csv")
    p.add_argument("--neighbor-summary-csv", required=True, help="stable_shell_neighbor_summary.csv")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def safe_mwu(a: pd.Series, b: pd.Series) -> tuple[float, float]:
    xa = pd.to_numeric(a, errors="coerce").dropna().to_numpy(dtype=float)
    xb = pd.to_numeric(b, errors="coerce").dropna().to_numpy(dtype=float)
    if len(xa) == 0 or len(xb) == 0:
        return np.nan, np.nan
    try:
        u, p = mannwhitneyu(xa, xb, alternative="two-sided")
        return float(u), float(p)
    except Exception:
        return np.nan, np.nan


def median_or_nan(x: pd.Series) -> float:
    z = pd.to_numeric(x, errors="coerce").dropna()
    return float(z.median()) if len(z) else np.nan


def resolve_axis_col(df: pd.DataFrame, axis: str) -> str:
    candidates = [f"{axis}_dom", f"{axis}_state", axis]
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"missing axis column for {axis}")


def mean_abs_corr_to_neighbors(sample_values: pd.Series, neighbor_df: pd.DataFrame) -> float:
    vals = []
    x = pd.to_numeric(sample_values, errors="coerce").to_numpy(dtype=float)
    for _, row in neighbor_df.iterrows():
        y = pd.to_numeric(row, errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() < 10:
            continue
        c = np.corrcoef(x[mask], y[mask])[0, 1]
        if np.isfinite(c):
            vals.append(abs(c))
    return float(np.mean(vals)) if vals else np.nan


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    merged = pd.read_csv(args.merged_curated_csv)
    merged["sample_id"] = merged["sample_id"].astype(str)

    mat = pd.read_parquet(args.selected_probe_matrix)
    if "ID_REF" not in mat.columns:
        raise SystemExit("selected probe matrix must contain ID_REF")
    probe_mat = mat.set_index("ID_REF")
    sample_cols = probe_mat.columns.astype(str).tolist()
    probe_mat.columns = sample_cols

    shell = pd.read_csv(args.stable_shell_csv)
    shell["sample_id"] = shell["sample_id"].astype(str)

    neighbor = pd.read_csv(args.neighbor_csv)
    neighbor_summary = pd.read_csv(args.neighbor_summary_csv)

    for col in ["sample_id", "state_domain"]:
        if col not in merged.columns:
            raise SystemExit(f"missing required column in merged curated table: {col}")

    d3_ids = merged.loc[merged["state_domain"] == "D3", "sample_id"].astype(str).tolist()
    non_d3_ids = merged.loc[merged["state_domain"] != "D3", "sample_id"].astype(str).tolist()

    d3_ids = [s for s in d3_ids if s in sample_cols]
    non_d3_ids = [s for s in non_d3_ids if s in sample_cols]

    sample_metrics = []
    shell_ids = shell["sample_id"].astype(str).tolist()

    axis_cols = {ax: resolve_axis_col(merged, ax) for ax in ["H", "S", "M", "R", "phi"]}

    for sid in [s for s in merged["sample_id"].astype(str).tolist() if s in sample_cols]:
        vals = pd.to_numeric(probe_mat[sid], errors="coerce").to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]

        mean_beta = float(np.mean(vals)) if len(vals) else np.nan
        std_beta = float(np.std(vals)) if len(vals) else np.nan
        iqr_beta = float(np.quantile(vals, 0.75) - np.quantile(vals, 0.25)) if len(vals) else np.nan
        mad_beta = float(np.median(np.abs(vals - np.median(vals)))) if len(vals) else np.nan

        row = merged.loc[merged["sample_id"] == sid].iloc[0]
        hrsm = np.array([row[axis_cols["H"]], row[axis_cols["S"]], row[axis_cols["M"]], row[axis_cols["R"]], row[axis_cols["phi"]]], dtype=float)

        min_shell_hrsm_dist = np.nan
        for sh in shell_ids:
            if sh == sid or sh not in merged["sample_id"].astype(str).tolist():
                continue
            sh_row = merged.loc[merged["sample_id"] == sh].iloc[0]
            sh_vec = np.array([sh_row[axis_cols["H"]], sh_row[axis_cols["S"]], sh_row[axis_cols["M"]], sh_row[axis_cols["R"]], sh_row[axis_cols["phi"]]], dtype=float)
            d = np.linalg.norm(hrsm - sh_vec)
            if not np.isfinite(min_shell_hrsm_dist) or d < min_shell_hrsm_dist:
                min_shell_hrsm_dist = float(d)

        if sid in probe_mat.columns:
            x = probe_mat[sid].to_numpy(dtype=float)
            min_shell_probe_dist = np.nan
            for sh in shell_ids:
                if sh == sid or sh not in probe_mat.columns:
                    continue
                y = probe_mat[sh].to_numpy(dtype=float)
                mask = np.isfinite(x) & np.isfinite(y)
                if mask.sum() == 0:
                    continue
                d = np.linalg.norm(x[mask] - y[mask]) / np.sqrt(mask.sum())
                if not np.isfinite(min_shell_probe_dist) or d < min_shell_probe_dist:
                    min_shell_probe_dist = float(d)
        else:
            min_shell_probe_dist = np.nan

        sample_metrics.append(
            {
                "sample_id": sid,
                "state_domain": row["state_domain"],
                "placeholder_condition_cur": row["placeholder_condition_cur"] if "placeholder_condition_cur" in row.index else np.nan,
                "tumor_status": row["tumor_status"] if "tumor_status" in row.index else np.nan,
                "candidate_stable_shell": bool(row["candidate_stable_shell"]) if "candidate_stable_shell" in row.index and pd.notna(row["candidate_stable_shell"]) else False,
                "mean_beta": mean_beta,
                "std_beta": std_beta,
                "iqr_beta": iqr_beta,
                "mad_beta": mad_beta,
                "mean_neighbor_distance": row["mean_neighbor_distance"] if "mean_neighbor_distance" in row.index else np.nan,
                "H": row[axis_cols["H"]],
                "S": row[axis_cols["S"]],
                "M": row[axis_cols["M"]],
                "R": row[axis_cols["R"]],
                "phi": row[axis_cols["phi"]],
                "min_hrsm_distance_to_shell": min_shell_hrsm_dist,
                "min_probe_distance_to_shell": min_shell_probe_dist,
            }
        )

    sample_metrics_df = pd.DataFrame(sample_metrics)

    tests = []
    for metric in [
        "mean_beta",
        "std_beta",
        "iqr_beta",
        "mad_beta",
        "mean_neighbor_distance",
        "min_hrsm_distance_to_shell",
        "min_probe_distance_to_shell",
        "H",
        "S",
        "M",
        "R",
        "phi",
    ]:
        a = sample_metrics_df.loc[sample_metrics_df["state_domain"] == "D3", metric]
        b = sample_metrics_df.loc[sample_metrics_df["state_domain"] != "D3", metric]
        u, p = safe_mwu(a, b)
        tests.append(
            {
                "comparison": "D3_vs_not_D3",
                "metric": metric,
                "n_D3": int(a.notna().sum()),
                "n_not_D3": int(b.notna().sum()),
                "median_D3": median_or_nan(a),
                "median_not_D3": median_or_nan(b),
                "delta_median_D3_minus_not_D3": median_or_nan(a) - median_or_nan(b),
                "u_statistic": u,
                "p_value": p,
            }
        )

    tests_df = pd.DataFrame(tests)

    curated = sample_metrics_df[sample_metrics_df["placeholder_condition_cur"].isin(["mgmt_methylated", "mgmt_non_methylated"])].copy()
    d3_meth = int(((curated["state_domain"] == "D3") & (curated["placeholder_condition_cur"] == "mgmt_methylated")).sum())
    d3_non = int(((curated["state_domain"] == "D3") & (curated["placeholder_condition_cur"] == "mgmt_non_methylated")).sum())
    other_meth = int(((curated["state_domain"] != "D3") & (curated["placeholder_condition_cur"] == "mgmt_methylated")).sum())
    other_non = int(((curated["state_domain"] != "D3") & (curated["placeholder_condition_cur"] == "mgmt_non_methylated")).sum())

    try:
        odds, p_fisher = fisher_exact([[d3_meth, d3_non], [other_meth, other_non]])
    except Exception:
        odds, p_fisher = np.nan, np.nan

    mgmt_df = pd.DataFrame(
        [
            {
                "focus_domain": "D3",
                "meth_in_domain": d3_meth,
                "nonmeth_in_domain": d3_non,
                "meth_outside_domain": other_meth,
                "nonmeth_outside_domain": other_non,
                "frac_meth_in_domain": d3_meth / max(d3_meth + d3_non, 1),
                "frac_meth_outside_domain": other_meth / max(other_meth + other_non, 1),
                "fisher_odds_ratio": odds,
                "fisher_p_value": p_fisher,
            }
        ]
    )

    # D3 sample proximity to stable shell
    d3_shell_df = sample_metrics_df.loc[sample_metrics_df["state_domain"] == "D3"].copy()
    d3_shell_df = d3_shell_df.sort_values(["candidate_stable_shell", "min_hrsm_distance_to_shell"], ascending=[False, True]).reset_index(drop=True)

    # Shell neighborhood overlap with D3
    if {"shell_sample_id", "neighbor_sample_id"}.issubset(neighbor.columns):
        d3_lookup = set(d3_ids)
        neighbor = neighbor.copy()
        neighbor["neighbor_sample_id"] = neighbor["neighbor_sample_id"].astype(str)
        neighbor["neighbor_is_D3"] = neighbor["neighbor_sample_id"].isin(d3_lookup)
        d3_overlap = (
            neighbor.groupby("shell_sample_id", as_index=False)["neighbor_is_D3"]
            .agg(["sum", "count"])
            .reset_index()
        )
        d3_overlap.columns = ["shell_sample_id", "n_D3_neighbors", "n_total_neighbors"]
        d3_overlap["frac_D3_neighbors"] = d3_overlap["n_D3_neighbors"] / d3_overlap["n_total_neighbors"]
    else:
        d3_overlap = pd.DataFrame(columns=["shell_sample_id", "n_D3_neighbors", "n_total_neighbors", "frac_D3_neighbors"])

    sample_metrics_df.to_csv(outdir / "d3_sample_level_metrics.csv", index=False)
    tests_df.to_csv(outdir / "d3_metric_tests.csv", index=False)
    mgmt_df.to_csv(outdir / "d3_mgmt_enrichment.csv", index=False)
    d3_shell_df.to_csv(outdir / "d3_proximity_to_stable_shell.csv", index=False)
    d3_overlap.to_csv(outdir / "d3_overlap_with_shell_neighborhoods.csv", index=False)

    summary = {
        "focus": "D3 audit: level vs stochasticity, MGMT enrichment, stable-shell proximity, and local neighborhood position",
        "n_total_samples": int(len(sample_metrics_df)),
        "n_D3": int((sample_metrics_df["state_domain"] == "D3").sum()),
        "n_not_D3": int((sample_metrics_df["state_domain"] != "D3").sum()),
        "n_D3_stable_shell": int(((sample_metrics_df["state_domain"] == "D3") & (sample_metrics_df["candidate_stable_shell"])).sum()),
    }
    with open(outdir / "d3_state_audit_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote D3 audit outputs to", outdir)
    print("[info] summary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
