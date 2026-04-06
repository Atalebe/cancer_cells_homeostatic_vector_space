#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import pandas as pd


def split_gene_tokens(val: str):
    if pd.isna(val):
        return []
    txt = str(val).strip()
    if txt == "":
        return []
    for sep in [";", ",", "///"]:
        txt = txt.replace(sep, "|")
    parts = [x.strip() for x in txt.split("|")]
    return [x for x in parts if x and x.upper() not in {"NA", "N/A", "NULL"}]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotated-csv", required=True)
    ap.add_argument("--outcsv", required=True)
    ap.add_argument("--outjson", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.annotated_csv, low_memory=False)

    gene_col = None
    for c in ["Gene_Symbol", "UCSC_RefGene_Name", "Gene", "Symbol"]:
        if c in df.columns:
            gene_col = c
            break

    counts = Counter()
    if gene_col is not None:
        for val in df[gene_col]:
            for tok in split_gene_tokens(val):
                counts[tok] += 1

    out = pd.DataFrame(
        [{"gene": k, "count": v} for k, v in counts.most_common(100)]
    )
    out.to_csv(args.outcsv, index=False)

    summary = {
        "annotated_csv": args.annotated_csv,
        "gene_col_used": gene_col,
        "n_rows": int(len(df)),
        "n_unique_gene_tokens": int(len(counts)),
        "top_20_gene_tokens": counts.most_common(20),
    }
    with open(args.outjson, "w") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote gene summary:", args.outcsv)
    print("[ok] wrote summary:", args.outjson)
    print("[info] gene_col_used:", gene_col)
    print("[info] n_unique_gene_tokens:", len(counts))


if __name__ == "__main__":
    main()
