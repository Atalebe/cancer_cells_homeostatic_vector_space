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
    "gene_name",
    "gene_names",
]

DIRECTION_COL_CANDIDATES = [
    "direction",
    "change_direction",
    "contrast_direction",
]

DELIM_RE = re.compile(r"[;|/,]")


def pick_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def normalize_gene_token(x: str) -> str:
    return x.strip().upper()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-csv", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--label", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input_csv)
    gene_col = pick_col(df, GENE_COL_CANDIDATES)
    direction_col = pick_col(df, DIRECTION_COL_CANDIDATES)

    token_rows = []
    if gene_col is not None:
        for _, row in df.iterrows():
            raw = row.get(gene_col)
            if pd.isna(raw):
                continue
            for tok in DELIM_RE.split(str(raw)):
                tok = normalize_gene_token(tok)
                if tok and tok not in {"NA", "NAN", "."}:
                    token_rows.append({
                        "gene_token": tok,
                        "direction": row.get(direction_col, pd.NA) if direction_col else pd.NA,
                        "ID_REF": row.get("ID_REF", pd.NA),
                    })

    long_df = pd.DataFrame(token_rows)

    if long_df.empty:
        counts = pd.DataFrame(columns=["gene_token", "n_probes"])
        by_direction = pd.DataFrame(columns=["gene_token", "direction", "n_probes"])
    else:
        counts = (
            long_df.groupby("gene_token")
            .size()
            .reset_index(name="n_probes")
            .sort_values(["n_probes", "gene_token"], ascending=[False, True])
            .reset_index(drop=True)
        )

        if direction_col is not None:
            by_direction = (
                long_df.groupby(["gene_token", "direction"], dropna=False)
                .size()
                .reset_index(name="n_probes")
                .sort_values(["n_probes", "gene_token"], ascending=[False, True])
                .reset_index(drop=True)
            )
        else:
            by_direction = pd.DataFrame(columns=["gene_token", "direction", "n_probes"])

    counts.to_csv(outdir / f"{args.label}_gene_token_counts.csv", index=False)
    by_direction.to_csv(outdir / f"{args.label}_gene_token_by_direction.csv", index=False)

    summary = {
        "input_csv": args.input_csv,
        "label": args.label,
        "gene_col_used": gene_col,
        "direction_col_used": direction_col,
        "n_rows_input": int(len(df)),
        "n_gene_tokens": int(len(long_df)),
        "n_unique_gene_tokens": int(counts["gene_token"].nunique()) if not counts.empty else 0,
        "top_20_gene_tokens": counts.head(20)["gene_token"].tolist() if not counts.empty else [],
    }

    with open(outdir / f"{args.label}_gene_symbol_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    if gene_col is None:
        print(f"[warn] no gene-like column found for {args.label}, wrote empty gene summaries")
    else:
        print(f"[ok] wrote gene summaries for {args.label}")
        print(f"[info] gene_col_used: {gene_col}")


if __name__ == "__main__":
    main()
