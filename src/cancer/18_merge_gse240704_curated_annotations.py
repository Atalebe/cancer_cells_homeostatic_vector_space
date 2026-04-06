from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge curated biological metadata into GSE240704 sample annotations."
    )
    parser.add_argument(
        "--sample-annotations",
        required=True,
        help="Existing sample_annotations parquet.",
    )
    parser.add_argument(
        "--curated-mapping",
        required=True,
        help="Curated metadata mapping CSV with sample_id and biology columns.",
    )
    parser.add_argument(
        "--outparquet",
        required=True,
        help="Output curated sample annotations parquet.",
    )
    parser.add_argument(
        "--outcsv",
        required=True,
        help="Output curated sample annotations CSV.",
    )
    return parser.parse_args()


def choose_first_present(row: pd.Series, candidates: list[str]):
    for c in candidates:
        if c in row.index and pd.notna(row[c]) and str(row[c]).strip() != "":
            return row[c]
    return pd.NA


def main() -> None:
    args = parse_args()

    ann = pd.read_parquet(args.sample_annotations)
    cur = pd.read_csv(args.curated_mapping, dtype=str)

    if "sample_id" not in ann.columns or "sample_id" not in cur.columns:
        raise ValueError("Both inputs must contain sample_id")

    merged = ann.merge(cur, on="sample_id", how="left", suffixes=("", "_cur"))

    merged["placeholder_condition"] = merged.apply(
        lambda r: choose_first_present(r, [
            "biological_condition", "condition", "disease_group",
            "tumor_status", "cell_type"
        ]),
        axis=1,
    )

    merged["placeholder_group"] = merged.apply(
        lambda r: choose_first_present(r, [
            "disease_group", "cell_type", "patient_id"
        ]),
        axis=1,
    )

    merged["placeholder_batch"] = merged.apply(
        lambda r: choose_first_present(r, [
            "batch", "timepoint"
        ]),
        axis=1,
    )

    merged["placeholder_note"] = merged.apply(
        lambda r: choose_first_present(r, [
            "source_note", "sex", "age"
        ]),
        axis=1,
    )

    Path(args.outparquet).parent.mkdir(parents=True, exist_ok=True)
    merged.to_parquet(args.outparquet, index=False)
    merged.to_csv(args.outcsv, index=False)

    print("[ok] wrote curated annotations parquet:", args.outparquet)
    print("[ok] wrote curated annotations csv:", args.outcsv)

    for col in [
        "placeholder_condition",
        "placeholder_group",
        "placeholder_batch",
        "placeholder_note",
    ]:
        non_null = int(merged[col].notna().sum())
        uniq = int(merged[col].dropna().astype(str).nunique())
        print(f"[info] {col}: non_null={non_null} unique={uniq}")


if __name__ == "__main__":
    main()
