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

REFGENE_CANDIDATES = [
    "refgene_group",
    "UCSC_RefGene_Group",
]

CPG_CANDIDATES = [
    "relation_to_cpg_island",
    "Relation_to_UCSC_CpG_Island",
]

DIRECTION_CANDIDATES = [
    "direction",
]

WEIGHT_CANDIDATES = [
    "abs_delta_median",
    "abs_rank_biserial",
    "q_value_bh",
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


def build_gene_table(
    df: pd.DataFrame,
    direction_col: str,
    gene_col: str,
    weight_col: str | None,
) -> pd.DataFrame:
    rows: list[dict] = []
    for _, row in df.iterrows():
        direction = row[direction_col]
        weight = float(row[weight_col]) if weight_col and pd.notna(row[weight_col]) else 1.0
        for token in split_gene_tokens(row[gene_col]):
            rows.append(
                {
                    "direction": direction,
                    "gene_token": token,
                    "weight": weight,
                }
            )
    if not rows:
        return pd.DataFrame(columns=["direction", "gene_token", "n", "weighted_sum"])
    out = (
        pd.DataFrame(rows)
        .groupby(["direction", "gene_token"], as_index=False)
        .agg(n=("gene_token", "size"), weighted_sum=("weight", "sum"))
        .sort_values(["direction", "weighted_sum", "n", "gene_token"], ascending=[True, False, False, True])
        .reset_index(drop=True)
    )
    return out


def build_feature_table(
    df: pd.DataFrame,
    direction_col: str,
    feature_col: str,
) -> pd.DataFrame:
    sub = df[[direction_col, feature_col]].copy()
    sub[feature_col] = sub[feature_col].fillna("").astype(str).str.strip()
    sub = sub[sub[feature_col] != ""]
    if sub.empty:
        return pd.DataFrame(columns=["direction", feature_col, "n", "fraction_within_direction"])
    counts = (
        sub.groupby([direction_col, feature_col], as_index=False)
        .size()
        .rename(columns={"size": "n"})
    )
    totals = counts.groupby(direction_col)["n"].sum().to_dict()
    counts["fraction_within_direction"] = counts.apply(
        lambda r: r["n"] / totals[r[direction_col]] if totals.get(r[direction_col], 0) > 0 else 0.0,
        axis=1,
    )
    counts = counts.sort_values([direction_col, "n", feature_col], ascending=[True, False, True]).reset_index(drop=True)
    return counts


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-csv", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--label", default="D3_vs_not_D3")
    ap.add_argument("--top-n", type=int, default=30)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input_csv)

    direction_col = choose_col(df, DIRECTION_CANDIDATES)
    gene_col = choose_col(df, GENE_CANDIDATES)
    refgene_col = choose_col(df, REFGENE_CANDIDATES)
    cpg_col = choose_col(df, CPG_CANDIDATES)
    weight_col = choose_col(df, WEIGHT_CANDIDATES)

    if direction_col is None:
        raise SystemExit("missing direction column")
    if gene_col is None:
        raise SystemExit("missing gene column")

    gene_table = build_gene_table(df, direction_col, gene_col, weight_col)
    gene_top = (
        gene_table.groupby("direction", group_keys=False)
        .head(args.top_n)
        .reset_index(drop=True)
    )
    gene_top.to_csv(outdir / f"{args.label}_directional_gene_top{args.top_n}.csv", index=False)

    if refgene_col is not None:
        refgene_table = build_feature_table(df, direction_col, refgene_col)
        refgene_table.to_csv(outdir / f"{args.label}_directional_refgene_group_summary.csv", index=False)
    else:
        refgene_table = pd.DataFrame()

    if cpg_col is not None:
        cpg_table = build_feature_table(df, direction_col, cpg_col)
        cpg_table.to_csv(outdir / f"{args.label}_directional_cpg_relation_summary.csv", index=False)
    else:
        cpg_table = pd.DataFrame()

    compact_rows: list[dict] = []
    directions = list(df[direction_col].dropna().astype(str).unique())
    for d in directions:
        sub_gene = gene_top[gene_top["direction"] == d].head(10)
        sub_ref = refgene_table[refgene_table["direction"] == d].head(5) if not refgene_table.empty else pd.DataFrame()
        sub_cpg = cpg_table[cpg_table["direction"] == d].head(5) if not cpg_table.empty else pd.DataFrame()
        compact_rows.append(
            {
                "direction": d,
                "n_probes": int((df[direction_col].astype(str) == d).sum()),
                "top_gene_tokens": "; ".join(sub_gene["gene_token"].astype(str).tolist()),
                "top_refgene_contexts": "; ".join(sub_ref.iloc[:, 1].astype(str).tolist()) if not sub_ref.empty else "",
                "top_cpg_contexts": "; ".join(sub_cpg.iloc[:, 1].astype(str).tolist()) if not sub_cpg.empty else "",
            }
        )

    compact_df = pd.DataFrame(compact_rows)
    compact_df.to_csv(outdir / f"{args.label}_directional_compact_report.csv", index=False)

    summary = {
        "input_csv": args.input_csv,
        "label": args.label,
        "direction_col_used": direction_col,
        "gene_col_used": gene_col,
        "refgene_col_used": refgene_col,
        "cpg_col_used": cpg_col,
        "weight_col_used": weight_col,
        "n_rows_input": int(len(df)),
        "directions_present": directions,
        "top_n": int(args.top_n),
    }
    with open(outdir / f"{args.label}_directional_compact_report_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"[ok] wrote compact D3 directional report to {outdir}")
    print("[info] summary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
