#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def safe_str(x):
    if pd.isna(x):
        return ""
    return str(x)


def detect_probe_like_columns(df: pd.DataFrame):
    out = []
    for col in df.columns:
        s = df[col].dropna().astype(str)
        if s.empty:
            continue
        head = s.head(200)
        frac_cg = head.str.match(r"^cg\d+$", na=False).mean()
        frac_rs = head.str.match(r"^rs\d+$", na=False).mean()
        frac_ilmn = head.str.contains("Ilmn|ILMN", regex=True, na=False).mean()
        if max(frac_cg, frac_rs, frac_ilmn) > 0.10:
            out.append(
                {
                    "column": col,
                    "frac_cg_like": float(frac_cg),
                    "frac_rs_like": float(frac_rs),
                    "frac_ilmn_like": float(frac_ilmn),
                }
            )
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ann_path = Path(args.annotation)
    if not ann_path.exists():
        raise SystemExit(f"missing input: {ann_path}")

    df = pd.read_parquet(ann_path)

    summary = {
        "annotation_file": str(ann_path),
        "n_rows": int(len(df)),
        "n_columns": int(len(df.columns)),
        "columns": list(map(str, df.columns)),
    }

    pd.DataFrame({"column": df.columns, "dtype": [str(df[c].dtype) for c in df.columns]}).to_csv(
        outdir / "annotation_columns.csv", index=False
    )

    probe_like = detect_probe_like_columns(df)
    pd.DataFrame(probe_like).to_csv(outdir / "annotation_probe_like_columns.csv", index=False)

    preview_cols = list(df.columns[: min(20, len(df.columns))])
    df.head(30)[preview_cols].to_csv(outdir / "annotation_preview_head30.csv", index=False)

    candidate_text_cols = []
    for col in df.columns:
        non_null = df[col].dropna()
        if non_null.empty:
            continue
        if df[col].dtype == "object":
            candidate_text_cols.append(
                {
                    "column": col,
                    "non_null": int(non_null.shape[0]),
                    "n_unique_non_null": int(non_null.astype(str).nunique()),
                    "example_values": json.dumps(list(map(safe_str, non_null.astype(str).head(10).tolist()))),
                }
            )
    pd.DataFrame(candidate_text_cols).to_csv(outdir / "annotation_text_columns_summary.csv", index=False)

    with open(outdir / "annotation_reference_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote annotation inspection outputs to", outdir)
    print("[info] rows:", len(df))
    print("[info] columns:", len(df.columns))
    print("[info] inspect:")
    print(" ", outdir / "annotation_columns.csv")
    print(" ", outdir / "annotation_probe_like_columns.csv")
    print(" ", outdir / "annotation_preview_head30.csv")
    print(" ", outdir / "annotation_text_columns_summary.csv")
    print(" ", outdir / "annotation_reference_summary.json")


if __name__ == "__main__":
    main()
