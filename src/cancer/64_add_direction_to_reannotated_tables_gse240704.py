#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def choose_delta_col(df: pd.DataFrame) -> str | None:
    candidates = [
        "delta_median_a_minus_b",
        "delta_mean_a_minus_b",
        "delta_median",
        "delta_mean",
        "loading",
    ]
    for c in candidates:
        if c in df.columns:
            return c
    return None


def assign_direction(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    out = df.copy()

    def f(x: float) -> str:
        if pd.isna(x):
            return "unknown"
        if x > 0:
            return "higher_in_a"
        if x < 0:
            return "lower_in_a"
        return "no_change"

    out["direction"] = out[value_col].apply(f)
    out["direction_source_col"] = value_col
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    batch = []

    for inp in args.inputs:
        inpath = Path(inp)
        df = pd.read_csv(inpath)

        value_col = choose_delta_col(df)
        if value_col is None:
            raise SystemExit(f"no delta/loading column found in {inpath}")

        out = assign_direction(df, value_col)

        outname = inpath.stem.replace("_normalized_annotation", "") + "_with_direction.csv"
        outpath = outdir / outname
        out.to_csv(outpath, index=False)

        summary = {
            "input_csv": str(inpath),
            "output_csv": str(outpath),
            "direction_source_col": value_col,
            "n_rows": int(len(out)),
            "direction_counts": out["direction"].value_counts(dropna=False).to_dict(),
        }
        summary_path = outdir / f"{inpath.stem.replace('_normalized_annotation', '')}_with_direction_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2))

        print(f"[ok] wrote: {outpath}")
        print(f"[ok] wrote: {summary_path}")

        batch.append(summary)

    batch_path = outdir / "direction_batch_summary.json"
    batch_path.write_text(json.dumps(batch, indent=2))
    print(f"[ok] wrote batch summary: {batch_path}")


if __name__ == "__main__":
    main()
