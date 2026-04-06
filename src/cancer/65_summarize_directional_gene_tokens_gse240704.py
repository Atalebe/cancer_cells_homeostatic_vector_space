#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


GENE_COL_CANDIDATES = [
    "ucsc_refgene_name",
    "UCSC_RefGene_Name",
    "gene_symbol",
    "gene_symbols",
    "gene",
]

DIRECTION_COL_CANDIDATES = [
    "direction",
]

WEIGHT_COL_CANDIDATES = [
    "abs_delta_median",
    "abs_rank_biserial",
    "abs_loading",
    "abs_delta_mean",
]


def choose_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def split_gene_tokens(value: object) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    parts = re.split(r"[;,/|]+", text)
    tokens = []
    for p in parts:
        tok = p.strip()
        if tok and tok.lower() != "nan":
            tokens.append(tok)
    return tokens


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--topn", type=int, default=50)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input)

    gene_col = choose_col(df, GENE_COL_CANDIDATES)
    direction_col = choose_col(df, DIRECTION_COL_CANDIDATES)
    weight_col = choose_col(df, WEIGHT_COL_CANDIDATES)

    if gene_col is None or direction_col is None:
        summary = {
            "input_csv": args.input,
            "label": args.label,
            "gene_col_used": gene_col,
            "direction_col_used": direction_col,
            "weight_col_used": weight_col,
            "n_rows_input": int(len(df)),
            "note": "Required gene or direction column missing, skipped cleanly."
        }
        (outdir / f"{args.label}_directional_gene_summary.json").write_text(json.dumps(summary, indent=2))
        print("[warn] missing gene or direction column, wrote summary only")
        return

    rows = []
    for _, r in df.iterrows():
        direction = r.get(direction_col, "unknown")
        weight = r.get(weight_col, 1.0) if weight_col else 1.0
        try:
            weight = float(weight)
        except Exception:
            weight = 1.0
        for tok in split_gene_tokens(r.get(gene_col)):
            rows.append({
                "direction": direction,
                "gene_token": tok,
                "weight": weight,
            })

    long_df = pd.DataFrame(rows)

    if long_df.empty:
        summary = {
            "input_csv": args.input,
            "label": args.label,
            "gene_col_used": gene_col,
            "direction_col_used": direction_col,
            "weight_col_used": weight_col,
            "n_rows_input": int(len(df)),
            "n_gene_tokens": 0,
            "n_unique_gene_tokens": 0,
        }
        (outdir / f"{args.label}_directional_gene_summary.json").write_text(json.dumps(summary, indent=2))
        print("[warn] no gene tokens found")
        return

    counts = (
        long_df.groupby(["direction", "gene_token"], as_index=False)
        .agg(
            n=("gene_token", "size"),
            weighted_sum=("weight", "sum"),
        )
        .sort_values(["direction", "weighted_sum", "n", "gene_token"], ascending=[True, False, False, True])
        .reset_index(drop=True)
    )

    counts.to_csv(outdir / f"{args.label}_directional_gene_token_counts.csv", index=False)

    top_rows = []
    for direction, sub in counts.groupby("direction"):
        for _, r in sub.head(args.topn).iterrows():
            top_rows.append(r.to_dict())
    pd.DataFrame(top_rows).to_csv(outdir / f"{args.label}_directional_gene_token_top.csv", index=False)

    summary = {
        "input_csv": args.input,
        "label": args.label,
        "gene_col_used": gene_col,
        "direction_col_used": direction_col,
        "weight_col_used": weight_col,
        "n_rows_input": int(len(df)),
        "n_gene_tokens": int(len(long_df)),
        "n_unique_gene_tokens": int(long_df["gene_token"].nunique()),
        "directions_present": sorted(long_df["direction"].dropna().unique().tolist()),
    }
    (outdir / f"{args.label}_directional_gene_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[ok] wrote directional gene summaries for {args.label}")


if __name__ == "__main__":
    main()
