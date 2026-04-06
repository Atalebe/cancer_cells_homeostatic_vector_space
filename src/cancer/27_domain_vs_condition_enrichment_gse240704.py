#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, kruskal


AXES = ["H", "S", "M", "R", "phi"]


def ensure_outdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def load_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.suffix == ".parquet":
        return pd.read_parquet(p)
    return pd.read_csv(p)


def choose_condition_column(df: pd.DataFrame) -> str:
    candidates = [
        "placeholder_condition_cur",
        "biological_condition",
        "tumor_status",
        "placeholder_condition",
        "placeholder_condition_ann",
    ]
    for col in candidates:
        if col in df.columns and df[col].notna().sum() > 0:
            return col
    raise ValueError(f"No usable condition column found. Available columns: {list(df.columns)}")


def normalize_condition_value(x: object) -> object:
    if pd.isna(x):
        return np.nan
    s = str(x).strip().lower()
    if s in {"mgmt_methylated", "methylated"}:
        return "mgmt_methylated"
    if s in {"mgmt_non_methylated", "non_methylated", "non-methylated"}:
        return "mgmt_non_methylated"
    if "non" in s and "methyl" in s:
        return "mgmt_non_methylated"
    if "methyl" in s:
        return "mgmt_methylated"
    return s


def fisher_or_nan(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    table = np.array([[a, b], [c, d]], dtype=int)
    if table.sum() == 0:
        return np.nan, np.nan
    try:
        odds, p = fisher_exact(table, alternative="two-sided")
        return float(odds), float(p)
    except Exception:
        return np.nan, np.nan


def median_or_nan(series: pd.Series) -> float:
    s = pd.to_numeric(series, errors="coerce").dropna()
    if len(s) == 0:
        return np.nan
    return float(s.median())


def mean_or_nan(series: pd.Series) -> float:
    s = pd.to_numeric(series, errors="coerce").dropna()
    if len(s) == 0:
        return np.nan
    return float(s.mean())


def resolve_axis_column(df: pd.DataFrame, axis: str) -> str:
    candidates = [
        axis,
        f"{axis}_dom",
        f"{axis}_ann",
        f"{axis}_state",
        f"{axis}_state_tbl",
    ]
    for col in candidates:
        if col in df.columns:
            return col
    raise KeyError(f"Could not resolve axis column for {axis}. Available columns: {list(df.columns)}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--state-domains-csv", required=True)
    ap.add_argument("--sample-annotations-curated", required=True)
    ap.add_argument("--stable-shell-csv", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    ensure_outdir(outdir)

    domains = pd.read_csv(args.state_domains_csv)
    ann = load_table(args.sample_annotations_curated)
    stable = pd.read_csv(args.stable_shell_csv)

    merged = domains.merge(ann, on="sample_id", how="left", suffixes=("_dom", "_ann"))

    axis_map = {axis: resolve_axis_column(merged, axis) for axis in AXES}
    condition_col = choose_condition_column(merged)
    merged["condition_clean"] = merged[condition_col].map(normalize_condition_value)

    merged["is_mgmt_methylated"] = merged["condition_clean"].eq("mgmt_methylated")
    merged["is_mgmt_non_methylated"] = merged["condition_clean"].eq("mgmt_non_methylated")
    merged["is_curated_mgmt"] = merged["condition_clean"].isin(
        ["mgmt_methylated", "mgmt_non_methylated"]
    )

    stable_ids = set(stable["sample_id"].astype(str).tolist())
    merged["is_stable_shell"] = merged["sample_id"].astype(str).isin(stable_ids)

    domain_counts = (
        merged.groupby("state_domain", dropna=False)
        .size()
        .reset_index(name="n_total")
        .sort_values(["state_domain"])
        .reset_index(drop=True)
    )

    curated_condition_counts = (
        merged.loc[merged["is_curated_mgmt"]]
        .groupby(["state_domain", "condition_clean"], dropna=False)
        .size()
        .reset_index(name="n")
        .sort_values(["state_domain", "condition_clean"])
        .reset_index(drop=True)
    )

    mgmt_subset = merged.loc[merged["is_curated_mgmt"]].copy()

    enrichment_rows: list[dict] = []
    for domain in sorted(mgmt_subset["state_domain"].dropna().unique()):
        in_domain = mgmt_subset["state_domain"].eq(domain)
        meth_in = int((in_domain & mgmt_subset["is_mgmt_methylated"]).sum())
        nonmeth_in = int((in_domain & mgmt_subset["is_mgmt_non_methylated"]).sum())
        meth_out = int((~in_domain & mgmt_subset["is_mgmt_methylated"]).sum())
        nonmeth_out = int((~in_domain & mgmt_subset["is_mgmt_non_methylated"]).sum())

        odds, p = fisher_or_nan(meth_in, nonmeth_in, meth_out, nonmeth_out)

        enrichment_rows.append(
            {
                "state_domain": domain,
                "meth_in_domain": meth_in,
                "nonmeth_in_domain": nonmeth_in,
                "meth_outside_domain": meth_out,
                "nonmeth_outside_domain": nonmeth_out,
                "frac_meth_in_domain": meth_in / (meth_in + nonmeth_in)
                if (meth_in + nonmeth_in) > 0
                else np.nan,
                "frac_meth_outside_domain": meth_out / (meth_out + nonmeth_out)
                if (meth_out + nonmeth_out) > 0
                else np.nan,
                "fisher_odds_ratio": odds,
                "fisher_p_value": p,
            }
        )

    domain_enrichment = pd.DataFrame(enrichment_rows)
    if not domain_enrichment.empty:
        domain_enrichment = domain_enrichment.sort_values("state_domain").reset_index(drop=True)

    domain_axis_rows: list[dict] = []
    for domain, sub in merged.groupby("state_domain", dropna=False):
        row = {"state_domain": domain, "n": int(len(sub))}
        for axis in AXES:
            col = axis_map[axis]
            row[f"{axis}_source_col"] = col
            row[f"{axis}_median"] = median_or_nan(sub[col])
            row[f"{axis}_mean"] = mean_or_nan(sub[col])
        domain_axis_rows.append(row)

    domain_axis_summary = pd.DataFrame(domain_axis_rows).sort_values("state_domain").reset_index(drop=True)

    kruskal_rows: list[dict] = []
    for axis in AXES:
        col = axis_map[axis]
        groups = []
        labels = []
        for domain, sub in merged.groupby("state_domain", dropna=False):
            vals = pd.to_numeric(sub[col], errors="coerce").dropna()
            if len(vals) > 0:
                groups.append(vals.to_numpy())
                labels.append(domain)
        if len(groups) >= 2:
            stat, p = kruskal(*groups)
            kruskal_rows.append(
                {
                    "axis": axis,
                    "source_col": col,
                    "n_domains_tested": len(groups),
                    "domains": json.dumps([str(x) for x in labels]),
                    "kruskal_statistic": float(stat),
                    "p_value": float(p),
                }
            )

    kruskal_df = pd.DataFrame(kruskal_rows)
    if not kruskal_df.empty:
        kruskal_df = kruskal_df.sort_values("axis").reset_index(drop=True)

    stable_rows: list[dict] = []
    curated = merged.loc[merged["is_curated_mgmt"]].copy()

    n_shell = int(curated["is_stable_shell"].sum())
    n_rest = int((~curated["is_stable_shell"]).sum())

    meth_shell = int((curated["is_stable_shell"] & curated["is_mgmt_methylated"]).sum())
    nonmeth_shell = int((curated["is_stable_shell"] & curated["is_mgmt_non_methylated"]).sum())
    meth_rest = int((~curated["is_stable_shell"] & curated["is_mgmt_methylated"]).sum())
    nonmeth_rest = int((~curated["is_stable_shell"] & curated["is_mgmt_non_methylated"]).sum())
    shell_odds, shell_p = fisher_or_nan(meth_shell, nonmeth_shell, meth_rest, nonmeth_rest)

    stable_rows.append(
        {
            "comparison": "stable_shell_vs_rest_mgmt_enrichment",
            "n_shell_curated": n_shell,
            "n_rest_curated": n_rest,
            "meth_shell": meth_shell,
            "nonmeth_shell": nonmeth_shell,
            "meth_rest": meth_rest,
            "nonmeth_rest": nonmeth_rest,
            "fisher_odds_ratio": shell_odds,
            "fisher_p_value": shell_p,
        }
    )

    for axis in AXES:
        col = axis_map[axis]
        a = pd.to_numeric(curated.loc[curated["is_stable_shell"], col], errors="coerce").dropna()
        b = pd.to_numeric(curated.loc[~curated["is_stable_shell"], col], errors="coerce").dropna()
        if len(a) == 0 or len(b) == 0:
            continue
        stat, p = mannwhitneyu(a, b, alternative="two-sided")
        stable_rows.append(
            {
                "comparison": f"stable_shell_vs_rest_axis_{axis}",
                "source_col": col,
                "n_shell_curated": int(len(a)),
                "n_rest_curated": int(len(b)),
                "shell_median": float(a.median()),
                "rest_median": float(b.median()),
                "delta_median_shell_minus_rest": float(a.median() - b.median()),
                "u_statistic": float(stat),
                "p_value": float(p),
            }
        )

    stable_vs_rest = pd.DataFrame(stable_rows)

    shell_domain_counts = (
        merged.loc[merged["is_stable_shell"]]
        .groupby("state_domain", dropna=False)
        .size()
        .reset_index(name="n_stable_shell")
        .sort_values("state_domain")
        .reset_index(drop=True)
    )

    domain_counts.to_csv(outdir / "domain_counts_from_merged_annotations.csv", index=False)
    curated_condition_counts.to_csv(outdir / "domain_condition_counts_curated.csv", index=False)
    domain_enrichment.to_csv(outdir / "domain_mgmt_enrichment_fisher.csv", index=False)
    domain_axis_summary.to_csv(outdir / "domain_axis_summary.csv", index=False)
    kruskal_df.to_csv(outdir / "domain_axis_kruskal_tests.csv", index=False)
    stable_vs_rest.to_csv(outdir / "stable_shell_vs_rest_curated_tests.csv", index=False)
    shell_domain_counts.to_csv(outdir / "stable_shell_domain_counts.csv", index=False)

    summary = {
        "n_total_merged": int(len(merged)),
        "condition_column_used": condition_col,
        "axis_column_map": axis_map,
        "n_curated_mgmt_subset": int(merged["is_curated_mgmt"].sum()),
        "n_domains": int(merged["state_domain"].nunique(dropna=True)),
        "n_stable_shell_total": int(merged["is_stable_shell"].sum()),
        "stable_shell_sample_ids": sorted(list(stable_ids)),
    }
    with open(outdir / "domain_condition_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"[ok] wrote domain-condition outputs to {outdir}")
    print("[info] summary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
