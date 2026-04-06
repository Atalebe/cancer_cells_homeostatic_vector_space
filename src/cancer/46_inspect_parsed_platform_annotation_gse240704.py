#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotation-csv", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    p = Path(args.annotation_csv)
    if not p.exists():
        raise SystemExit(f"missing input: {p}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(p, dtype=str, low_memory=False)

    pd.DataFrame(
        {
            "column": df.columns,
            "non_null": [int(df[c].notna().sum()) for c in df.columns],
            "n_unique_non_null": [int(df[c].dropna().nunique()) for c in df.columns],
        }
    ).to_csv(outdir / "parsed_annotation_columns_summary.csv", index=False)

    df.head(40).to_csv(outdir / "parsed_annotation_preview_head40.csv", index=False)

    probe_like = []
    for col in df.columns:
        s = df[col].dropna().astype(str)
        if s.empty:
            continue
        frac_cg = s.head(500).str.match(r"^cg\d+$", na=False).mean()
        if frac_cg > 0.1:
            probe_like.append({"column": col, "frac_cg_like": float(frac_cg)})
    pd.DataFrame(probe_like).to_csv(outdir / "parsed_annotation_probe_like_columns.csv", index=False)

    summary = {
        "annotation_csv": str(p),
        "n_rows": int(len(df)),
        "n_columns": int(len(df.columns)),
        "columns": list(map(str, df.columns)),
    }
    with open(outdir / "parsed_annotation_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote parsed annotation inspection to", outdir)
    print("[info] inspect:")
    print(" ", outdir / "parsed_annotation_summary.json")
    print(" ", outdir / "parsed_annotation_columns_summary.csv")
    print(" ", outdir / "parsed_annotation_probe_like_columns.csv")
    print(" ", outdir / "parsed_annotation_preview_head40.csv")


if __name__ == "__main__":
    main()
