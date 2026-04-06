#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


DELIM_RE = re.compile(r"[;|/,]")


def extract_tokens(series: pd.Series) -> pd.DataFrame:
    rows = []
    for val in series.dropna():
        text = str(val).strip()
        if not text:
            continue
        for tok in DELIM_RE.split(text):
            tok = tok.strip().upper()
            if tok and tok not in {"NA", "NAN", "."}:
                rows.append({"gene_token": tok})
    return pd.DataFrame(rows)


def summarize_value_counts(df: pd.DataFrame, col: str) -> pd.DataFrame:
    s = df[col].astype("string").fillna("").str.strip()
    s = s[s != ""]
    if s.empty:
        return pd.DataFrame(columns=[col, "n", "fraction"])
    out = s.value_counts().rename_axis(col).reset_index(name="n")
    out["fraction"] = out["n"] / out["n"].sum()
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-csv", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input_csv)

    gene_col = None
    for c in ["ucsc_refgene_name", "UCSC_RefGene_Name", "gene_symbol", "gene_name"]:
        if c in df.columns:
            gene_col = c
            break

    if gene_col is not None:
        gene_tokens = extract_tokens(df[gene_col])
        if gene_tokens.empty:
            gene_counts = pd.DataFrame(columns=["gene_token", "n"])
        else:
            gene_counts = (
                gene_tokens.groupby("gene_token")
                .size()
                .reset_index(name="n")
                .sort_values(["n", "gene_token"], ascending=[False, True])
                .reset_index(drop=True)
            )
    else:
        gene_counts = pd.DataFrame(columns=["gene_token", "n"])

    gene_counts.to_csv(outdir / "s_axis_driver_gene_token_counts.csv", index=False)

    for col in [
        "refgene_group",
        "cpg_island_name",
        "relation_to_cpg_island",
        "regulatory_feature_group",
        "regulatory_feature_name",
        "phantom4_enhancers",
        "hmm_island",
    ]:
        if col in df.columns:
            summarize_value_counts(df, col).to_csv(outdir / f"s_axis_driver_{col}_summary.csv", index=False)
        else:
            pd.DataFrame(columns=[col, "n", "fraction"]).to_csv(
                outdir / f"s_axis_driver_{col}_summary.csv", index=False
            )

    summary = {
        "input_csv": args.input_csv,
        "n_rows_input": int(len(df)),
        "gene_col_used": gene_col,
        "n_unique_gene_tokens": int(gene_counts["gene_token"].nunique()) if not gene_counts.empty else 0,
        "top_20_gene_tokens": gene_counts.head(20)["gene_token"].tolist() if not gene_counts.empty else [],
    }

    with open(outdir / "s_axis_driver_biology_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote S-axis driver biology summaries")
    print(f"[info] gene_col_used: {gene_col}")
    print(f"[info] n_unique_gene_tokens: {summary['n_unique_gene_tokens']}")


if __name__ == "__main__":
    main()
