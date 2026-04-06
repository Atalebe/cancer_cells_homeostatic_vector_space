#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--outcsv", required=True)
    ap.add_argument("--summary-json", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.input)

    if "loading" not in df.columns:
        raise SystemExit("input file missing loading column")

    out = df.copy()

    def f(x: float) -> str:
        if pd.isna(x):
            return "unknown"
        if x > 0:
            return "positive_loading"
        if x < 0:
            return "negative_loading"
        return "zero_loading"

    out["direction"] = out["loading"].apply(f)
    out["direction_source_col"] = "loading"
    out.to_csv(args.outcsv, index=False)

    summary = {
        "input_csv": args.input,
        "output_csv": args.outcsv,
        "direction_source_col": "loading",
        "n_rows": int(len(out)),
        "direction_counts": out["direction"].value_counts(dropna=False).to_dict(),
    }
    Path(args.summary_json).write_text(json.dumps(summary, indent=2))

    print(f"[ok] wrote signed S-axis table: {args.outcsv}")
    print(f"[ok] wrote summary: {args.summary_json}")


if __name__ == "__main__":
    main()
