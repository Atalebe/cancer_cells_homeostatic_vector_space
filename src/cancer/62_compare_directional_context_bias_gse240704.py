#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


DIRECTION_COLS = ["direction", "change_direction", "contrast_direction"]
FEATURE_COLS = [
    "refgene_group",
    "cpg_island_name",
    "regulatory_feature_group",
    "regulatory_feature_name",
    "phantom4_enhancers",
    "hmm_island",
]


def pick_direction_col(df: pd.DataFrame) -> str | None:
    for c in DIRECTION_COLS:
        if c in df.columns:
            return c
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-csv", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--label", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input_csv)
    dcol = pick_direction_col(df)

    if dcol is None:
        summary = {
            "input_csv": args.input_csv,
            "label": args.label,
            "direction_col_used": None,
            "n_rows_input": int(len(df)),
            "note": "No direction column present, directional context bias skipped.",
        }
        with open(outdir / f"{args.label}_directional_context_bias_summary.json", "w") as f:
            json.dump(summary, f, indent=2)
        print(f"[warn] no direction column found for {args.label}, skipped directional context summaries")
        return

    summary = {
        "input_csv": args.input_csv,
        "label": args.label,
        "direction_col_used": dcol,
        "n_rows_input": int(len(df)),
        "feature_outputs": {},
    }

    for feat in FEATURE_COLS:
        if feat not in df.columns:
            pd.DataFrame(columns=[dcol, feat, "n", "fraction_within_direction"]).to_csv(
                outdir / f"{args.label}_{feat}_by_direction.csv", index=False
            )
            summary["feature_outputs"][feat] = {"present": False}
            continue

        sub = df[[dcol, feat]].copy()
        sub[dcol] = sub[dcol].astype("string").fillna("").str.strip()
        sub[feat] = sub[feat].astype("string").fillna("").str.strip()
        sub = sub[(sub[dcol] != "") & (sub[feat] != "")]

        if sub.empty:
            out = pd.DataFrame(columns=[dcol, feat, "n", "fraction_within_direction"])
        else:
            out = (
                sub.groupby([dcol, feat], dropna=False)
                .size()
                .reset_index(name="n")
            )
            totals = out.groupby(dcol)["n"].transform("sum")
            out["fraction_within_direction"] = out["n"] / totals
            out = out.sort_values([dcol, "n"], ascending=[True, False]).reset_index(drop=True)

        out.to_csv(outdir / f"{args.label}_{feat}_by_direction.csv", index=False)
        summary["feature_outputs"][feat] = {
            "present": True,
            "n_non_empty_rows": int(len(sub)),
            "n_unique_values": int(sub[feat].nunique()) if not sub.empty else 0,
        }

    with open(outdir / f"{args.label}_directional_context_bias_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[ok] wrote directional context bias summaries for {args.label}")
    print(f"[info] direction_col_used: {dcol}")


if __name__ == "__main__":
    main()
