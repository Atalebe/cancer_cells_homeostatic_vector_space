#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


AXES = ["H", "S", "M", "R", "phi"]


def median_or_nan(x: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce").dropna()
    if len(x) == 0:
        return float("nan")
    return float(np.median(x))


def mw(a: pd.Series, b: pd.Series):
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) == 0 or len(b) == 0:
        return np.nan, np.nan
    stat, p = mannwhitneyu(a, b, alternative="two-sided")
    return float(stat), float(p)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--merged-curated-csv", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.merged_curated_csv)

    rows = []
    for axis in AXES:
        col = f"{axis}_dom"
        for a, b in [("D3", "D1"), ("D3", "D2"), ("D2", "D1")]:
            xa = df.loc[df["state_domain"] == a, col]
            xb = df.loc[df["state_domain"] == b, col]
            u, p = mw(xa, xb)
            rows.append(
                {
                    "domain_a": a,
                    "domain_b": b,
                    "axis": axis,
                    "n_a": int(pd.to_numeric(xa, errors="coerce").notna().sum()),
                    "n_b": int(pd.to_numeric(xb, errors="coerce").notna().sum()),
                    "median_a": median_or_nan(xa),
                    "median_b": median_or_nan(xb),
                    "delta_median_a_minus_b": median_or_nan(xa) - median_or_nan(xb),
                    "u_statistic": u,
                    "p_value": p,
                }
            )

    out = pd.DataFrame(rows)
    out.to_csv(outdir / "domain_extreme_signature_tests.csv", index=False)

    summary = {
        "domains_present": sorted(df["state_domain"].dropna().astype(str).unique().tolist()),
        "n_rows": int(len(df)),
        "focus": "D3_vs_D1_D2_extreme_state_tests",
    }
    with open(outdir / "domain_extreme_signature_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"[ok] wrote extreme signature tests: {outdir / 'domain_extreme_signature_tests.csv'}")
    print(f"[ok] wrote summary: {outdir / 'domain_extreme_signature_summary.json'}")


if __name__ == "__main__":
    main()
