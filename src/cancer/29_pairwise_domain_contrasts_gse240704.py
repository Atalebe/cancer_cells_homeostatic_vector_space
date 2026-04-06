#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scipy.stats import mannwhitneyu
except Exception as exc:
    raise SystemExit(f"scipy is required: {exc}")


AXES = ["H", "S", "M", "R", "phi"]


def median_or_nan(x: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce").dropna()
    if len(x) == 0:
        return float("nan")
    return float(np.median(x))


def mean_or_nan(x: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce").dropna()
    if len(x) == 0:
        return float("nan")
    return float(np.mean(x))


def safe_mwu(a: pd.Series, b: pd.Series):
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) == 0 or len(b) == 0:
        return np.nan, np.nan, np.nan
    res = mannwhitneyu(a, b, alternative="two-sided")
    u = float(res.statistic)
    p = float(res.pvalue)
    rbc = (2.0 * u) / (len(a) * len(b)) - 1.0
    return u, p, float(rbc)


def resolve_axis_cols(df: pd.DataFrame) -> dict[str, str]:
    out = {}
    for axis in AXES:
        candidates = [axis, f"{axis}_dom", f"{axis}_state"]
        chosen = next((c for c in candidates if c in df.columns), None)
        if chosen is None:
            raise KeyError(f"could not resolve source column for axis '{axis}'")
        out[axis] = chosen
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain-followup-dir", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    in_dir = Path(args.domain_followup_dir)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    axis_summary_path = in_dir / "domain_axis_summary.csv"
    curated_overlay_path = in_dir / "domain_condition_merged_curated.csv"

    if not axis_summary_path.exists():
        raise SystemExit(f"missing input: {axis_summary_path}")
    if not curated_overlay_path.exists():
        raise SystemExit(f"missing input: {curated_overlay_path}")

    axis_summary = pd.read_csv(axis_summary_path)
    df = pd.read_csv(curated_overlay_path)

    axis_map = resolve_axis_cols(df)

    if "state_domain" not in df.columns:
        raise SystemExit("state_domain column not found in merged curated table")

    domain_counts = (
        df["state_domain"]
        .astype(str)
        .value_counts(dropna=False)
        .rename_axis("state_domain")
        .reset_index(name="n")
        .sort_values(["state_domain"])
        .reset_index(drop=True)
    )
    domain_counts.to_csv(out_dir / "pairwise_domain_counts.csv", index=False)

    domains = sorted([d for d in df["state_domain"].dropna().astype(str).unique() if d != "nan"])

    rows = []
    for d1, d2 in combinations(domains, 2):
        sub1 = df.loc[df["state_domain"].astype(str) == d1].copy()
        sub2 = df.loc[df["state_domain"].astype(str) == d2].copy()

        for axis in AXES:
            col = axis_map[axis]
            u, p, rbc = safe_mwu(sub1[col], sub2[col])

            rows.append(
                {
                    "domain_a": d1,
                    "domain_b": d2,
                    "axis": axis,
                    "source_col": col,
                    "n_a": int(pd.to_numeric(sub1[col], errors="coerce").notna().sum()),
                    "n_b": int(pd.to_numeric(sub2[col], errors="coerce").notna().sum()),
                    "median_a": median_or_nan(sub1[col]),
                    "median_b": median_or_nan(sub2[col]),
                    "mean_a": mean_or_nan(sub1[col]),
                    "mean_b": mean_or_nan(sub2[col]),
                    "delta_median_a_minus_b": median_or_nan(sub1[col]) - median_or_nan(sub2[col]),
                    "u_statistic": u,
                    "p_value": p,
                    "rank_biserial": rbc,
                }
            )

    out = pd.DataFrame(rows).sort_values(["axis", "domain_a", "domain_b"]).reset_index(drop=True)
    out.to_csv(out_dir / "pairwise_domain_contrasts.csv", index=False)

    summary = {
        "n_domains": len(domains),
        "domains": domains,
        "n_pairwise_tests": int(len(out)),
        "axis_column_map": axis_map,
    }
    with open(out_dir / "pairwise_domain_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"[ok] wrote pairwise domain contrasts to {out_dir}")
    print("[info] summary:")
    print(f"  n_domains: {len(domains)}")
    print(f"  domains: {domains}")
    print(f"  n_pairwise_tests: {len(out)}")


if __name__ == "__main__":
    main()
