#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def safe_mode(series: pd.Series) -> str:
    vals = series.dropna().astype(str)
    if vals.empty:
        return ""
    vc = vals.value_counts()
    return str(vc.index[0])


def first_matching_col(columns, candidates):
    for c in candidates:
        if c in columns:
            return c
    low_map = {str(x).lower(): x for x in columns}
    for c in candidates:
        if c.lower() in low_map:
            return low_map[c.lower()]
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotated-contrast-csv", required=True)
    ap.add_argument("--outcsv", required=True)
    ap.add_argument("--outjson", required=True)
    ap.add_argument("--top-n", type=int, default=200)
    args = ap.parse_args()

    df = pd.read_csv(args.annotated_contrast_csv)
    df = df.sort_values(["q_value_bh", "abs_delta_median", "abs_rank_biserial"], ascending=[True, False, False]).head(args.top_n).copy()

    gene_col = first_matching_col(df.columns, [
        "UCSC_RefGene_Name", "Gene_Name", "gene", "GeneSymbol", "genes"
    ])
    feature_col = first_matching_col(df.columns, [
        "UCSC_RefGene_Group", "Feature_Type", "feature", "GeneFeature"
    ])
    island_col = first_matching_col(df.columns, [
        "Relation_to_UCSC_CpG_Island", "CpG_Island", "island", "Island"
    ])
    chr_col = first_matching_col(df.columns, [
        "CHR", "chr", "Chromosome"
    ])

    rows = []

    for direction, sub in df.groupby("direction", dropna=False):
        row = {
            "direction": direction,
            "n_probes": int(len(sub)),
            "median_abs_delta": float(sub["abs_delta_median"].median()),
            "max_abs_delta": float(sub["abs_delta_median"].max()),
            "median_abs_rank_biserial": float(sub["abs_rank_biserial"].median()),
            "best_q_value_bh": float(sub["q_value_bh"].min()),
            "top_probe_example": str(sub.iloc[0]["ID_REF"]),
        }
        if gene_col:
            genes = (
                sub[gene_col]
                .dropna()
                .astype(str)
                .str.split(";|,")
                .explode()
                .str.strip()
            )
            genes = genes[genes != ""]
            row["top_gene_tokens_json"] = json.dumps(genes.value_counts().head(15).to_dict())
        else:
            row["top_gene_tokens_json"] = json.dumps({})

        if feature_col:
            row["dominant_feature"] = safe_mode(sub[feature_col])
            row["feature_counts_json"] = json.dumps(sub[feature_col].dropna().astype(str).value_counts().head(10).to_dict())
        else:
            row["dominant_feature"] = ""
            row["feature_counts_json"] = json.dumps({})

        if island_col:
            row["dominant_island_context"] = safe_mode(sub[island_col])
            row["island_counts_json"] = json.dumps(sub[island_col].dropna().astype(str).value_counts().head(10).to_dict())
        else:
            row["dominant_island_context"] = ""
            row["island_counts_json"] = json.dumps({})

        if chr_col:
            row["top_chromosomes_json"] = json.dumps(sub[chr_col].dropna().astype(str).value_counts().head(10).to_dict())
        else:
            row["top_chromosomes_json"] = json.dumps({})

        rows.append(row)

    out = pd.DataFrame(rows).sort_values("direction").reset_index(drop=True)
    Path(args.outcsv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.outcsv, index=False)

    summary = {
        "input_rows": int(len(pd.read_csv(args.annotated_contrast_csv))),
        "rows_used": int(len(df)),
        "top_n": int(args.top_n),
        "gene_col": gene_col,
        "feature_col": feature_col,
        "island_col": island_col,
        "chr_col": chr_col,
        "directions_present": sorted(out["direction"].astype(str).tolist()),
    }
    with open(args.outjson, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[ok] wrote signature summary: {args.outcsv}")
    print(f"[ok] wrote summary json: {args.outjson}")
    print("[info] summary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
