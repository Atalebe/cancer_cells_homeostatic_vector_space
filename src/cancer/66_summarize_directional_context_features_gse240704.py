#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


FEATURE_CANDIDATES = {
    "refgene_group": ["refgene_group", "UCSC_RefGene_Group"],
    "relation_to_cpg_island": ["relation_to_cpg_island", "Relation_to_UCSC_CpG_Island"],
    "regulatory_feature_group": ["regulatory_feature_group", "Regulatory_Feature_Group"],
    "regulatory_feature_name": ["regulatory_feature_name", "Regulatory_Feature_Name"],
    "hmm_island": ["hmm_island", "HMM_Island"],
}

DIRECTION_COL_CANDIDATES = ["direction"]


def choose_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def clean_series(s: pd.Series) -> pd.Series:
    x = s.astype(str).str.strip()
    x = x[(x != "") & (x.str.lower() != "nan")]
    return x


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input)
    direction_col = choose_col(df, DIRECTION_COL_CANDIDATES)

    if direction_col is None:
        summary = {
            "input_csv": args.input,
            "label": args.label,
            "direction_col_used": None,
            "note": "No direction column present, skipped cleanly."
        }
        (outdir / f"{args.label}_directional_context_summary.json").write_text(json.dumps(summary, indent=2))
        print("[warn] no direction column found")
        return

    feature_outputs = {}

    for feature_name, candidates in FEATURE_CANDIDATES.items():
        col = choose_col(df, candidates)
        if col is None:
            feature_outputs[feature_name] = {"column_used": None, "n_non_empty": 0}
            continue

        sub = df[[direction_col, col]].copy()
        sub[col] = clean_series(sub[col])
        sub = sub[sub[col].notna()]

        if sub.empty:
            feature_outputs[feature_name] = {"column_used": col, "n_non_empty": 0}
            continue

        counts = (
            sub.groupby([direction_col, col], as_index=False)
            .size()
            .rename(columns={"size": "n"})
        )
        totals = counts.groupby(direction_col)["n"].transform("sum")
        counts["fraction_within_direction"] = counts["n"] / totals
        counts = counts.sort_values([direction_col, "n", col], ascending=[True, False, True])

        out_csv = outdir / f"{args.label}_{feature_name}_by_direction.csv"
        counts.to_csv(out_csv, index=False)

        feature_outputs[feature_name] = {
            "column_used": col,
            "n_non_empty": int(len(sub)),
            "n_unique_non_empty": int(sub[col].nunique()),
            "output_csv": str(out_csv),
        }

    summary = {
        "input_csv": args.input,
        "label": args.label,
        "direction_col_used": direction_col,
        "n_rows_input": int(len(df)),
        "feature_outputs": feature_outputs,
    }
    (outdir / f"{args.label}_directional_context_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[ok] wrote directional context summaries for {args.label}")


if __name__ == "__main__":
    main()
