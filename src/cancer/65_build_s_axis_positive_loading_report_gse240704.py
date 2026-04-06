#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


GENE_CANDIDATES = [
    "ucsc_refgene_name",
    "UCSC_RefGene_Name",
    "gene_symbol",
    "Gene_Symbol",
]

CPG_CANDIDATES = [
    "relation_to_cpg_island",
    "Relation_to_UCSC_CpG_Island",
]

REFGENE_CANDIDATES = [
    "refgene_group",
    "UCSC_RefGene_Group",
]

DIRECTION_CANDIDATES = [
    "direction",
]

WEIGHT_CANDIDATES = [
    "abs_loading",
    "abs_delta_median",
]


def choose_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def split_gene_tokens(val: object) -> list[str]:
    if pd.isna(val):
        return []
    s = str(val).strip()
    if not s:
        return []
    tokens: list[str] = []
    for chunk in s.replace(",", ";").split(";"):
        token = chunk.strip()
        if token:
            tokens.append(token)
    return tokens


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-csv", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--positive-label", default="positive_loading")
    ap.add_argument("--top-n", type=int, default=40)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input_csv)

    direction_col = choose_col(df, DIRECTION_CANDIDATES)
    gene_col = choose_col(df, GENE_CANDIDATES)
    cpg_col = choose_col(df, CPG_CANDIDATES)
    refgene_col = choose_col(df, REFGENE_CANDIDATES)
    weight_col = choose_col(df, WEIGHT_CANDIDATES)

    if direction_col is None:
        raise SystemExit("missing direction column")
    if gene_col is None:
        raise SystemExit("missing gene column")

    pos = df[df[direction_col].astype(str) == args.positive_label].copy()
    if pos.empty:
        raise SystemExit(f"no rows found for direction={args.positive_label}")

    gene_rows: list[dict] = []
    for _, row in pos.iterrows():
        weight = float(row[weight_col]) if weight_col and pd.notna(row[weight_col]) else 1.0
        for token in split_gene_tokens(row[gene_col]):
            gene_rows.append({"gene_token": token, "weight": weight})

    if gene_rows:
        gene_df = (
            pd.DataFrame(gene_rows)
            .groupby("gene_token", as_index=False)
            .agg(n=("gene_token", "size"), weighted_sum=("weight", "sum"))
            .sort_values(["weighted_sum", "n", "gene_token"], ascending=[False, False, True])
            .reset_index(drop=True)
        )
    else:
        gene_df = pd.DataFrame(columns=["gene_token", "n", "weighted_sum"])

    gene_df.head(args.top_n).to_csv(outdir / "s_axis_positive_loading_top_genes.csv", index=False)

    top_probes = pos.copy()
    sort_col = weight_col if weight_col else pos.columns[0]
    top_probes = top_probes.sort_values(sort_col, ascending=False).head(args.top_n).reset_index(drop=True)
    top_probes.to_csv(outdir / "s_axis_positive_loading_top_probes.csv", index=False)

    if cpg_col and cpg_col in pos.columns:
        cpg_counts = (
            pos[cpg_col].fillna("").astype(str).str.strip().replace("", pd.NA).dropna()
            .value_counts(dropna=False)
            .rename_axis("relation_to_cpg_island")
            .reset_index(name="n")
        )
        total = int(cpg_counts["n"].sum()) if not cpg_counts.empty else 0
        if total > 0:
            cpg_counts["fraction"] = cpg_counts["n"] / total
        cpg_counts.to_csv(outdir / "s_axis_positive_loading_cpg_context_summary.csv", index=False)
    else:
        cpg_counts = pd.DataFrame()

    if refgene_col and refgene_col in pos.columns:
        ref_counts = (
            pos[refgene_col].fillna("").astype(str).str.strip().replace("", pd.NA).dropna()
            .value_counts(dropna=False)
            .rename_axis("refgene_group")
            .reset_index(name="n")
        )
        total = int(ref_counts["n"].sum()) if not ref_counts.empty else 0
        if total > 0:
            ref_counts["fraction"] = ref_counts["n"] / total
        ref_counts.to_csv(outdir / "s_axis_positive_loading_refgene_summary.csv", index=False)
    else:
        ref_counts = pd.DataFrame()

    summary = {
        "input_csv": args.input_csv,
        "direction_col_used": direction_col,
        "gene_col_used": gene_col,
        "cpg_col_used": cpg_col,
        "refgene_col_used": refgene_col,
        "weight_col_used": weight_col,
        "positive_label": args.positive_label,
        "n_rows_input": int(len(df)),
        "n_positive_rows": int(len(pos)),
        "n_unique_gene_tokens_positive": int(gene_df["gene_token"].nunique()) if not gene_df.empty else 0,
        "top_n": int(args.top_n),
    }
    with open(outdir / "s_axis_positive_loading_report_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"[ok] wrote S-axis positive loading report to {outdir}")
    print("[info] summary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
