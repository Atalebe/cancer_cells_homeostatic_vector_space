#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


FEATURE_GROUPS = {
    "refgene_group": ["refgene_group", "UCSC_RefGene_Group"],
    "ucsc_refgene_name": ["ucsc_refgene_name", "UCSC_RefGene_Name"],
    "ucsc_refgene_accession": ["ucsc_refgene_accession", "UCSC_RefGene_Accession"],
    "cpg_island_name": ["cpg_island_name", "UCSC_CpG_Islands_Name"],
    "relation_to_cpg_island": ["relation_to_cpg_island", "Relation_to_UCSC_CpG_Island"],
    "regulatory_feature_group": ["regulatory_feature_group", "Regulatory_Feature_Group"],
    "regulatory_feature_name": ["regulatory_feature_name", "Regulatory_Feature_Name"],
    "phantom4_enhancers": ["phantom4_enhancers", "Phantom4_Enhancers"],
    "hmm_island": ["hmm_island", "HMM_Island"],
    "chromosome": ["chromosome", "CHR"],
}


def pick_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def summarize_column(df: pd.DataFrame, col: str) -> pd.DataFrame:
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
    ap.add_argument("--label", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input_csv)

    summary = {
        "input_csv": args.input_csv,
        "label": args.label,
        "n_rows_input": int(len(df)),
        "columns_present": list(df.columns),
        "feature_outputs": {},
    }

    for group_name, candidates in FEATURE_GROUPS.items():
        col = pick_col(df, candidates)
        outpath = outdir / f"{args.label}_{group_name}_summary.csv"

        if col is None:
            pd.DataFrame(columns=[group_name, "n", "fraction"]).to_csv(outpath, index=False)
            summary["feature_outputs"][group_name] = {"column_used": None, "n_non_empty": 0}
            continue

        out = summarize_column(df, col)
        out.to_csv(outpath, index=False)

        summary["feature_outputs"][group_name] = {
            "column_used": col,
            "n_non_empty": int(out["n"].sum()) if not out.empty else 0,
            "n_unique_non_empty": int(out[col].nunique()) if not out.empty else 0,
        }

    with open(outdir / f"{args.label}_context_feature_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[ok] wrote context summaries for {args.label}")


if __name__ == "__main__":
    main()
