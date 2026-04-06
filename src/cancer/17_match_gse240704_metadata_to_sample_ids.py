from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Match extracted raw metadata rows to canonical sample_id values."
    )
    parser.add_argument(
        "--sample-registry",
        required=True,
        help="Path to sample registry CSV.",
    )
    parser.add_argument(
        "--metadata-csv",
        required=True,
        help="Path to extracted metadata candidate CSV.",
    )
    parser.add_argument(
        "--outcsv",
        required=True,
        help="Output matched metadata CSV.",
    )
    return parser.parse_args()


def sample_number_from_text(x: str) -> int | None:
    if pd.isna(x):
        return None
    text = str(x).strip()
    m = re.search(r"(?:sample)\s*([0-9]+)", text, flags=re.I)
    if m:
        return int(m.group(1))
    m = re.fullmatch(r"([0-9]+)", text)
    if m:
        return int(m.group(1))
    return None


def main() -> None:
    args = parse_args()

    reg = pd.read_csv(args.sample_registry)
    if "sample_id" not in reg.columns or "sample_number" not in reg.columns:
        raise ValueError("sample_registry must have sample_id and sample_number")

    meta = pd.read_csv(args.metadata_csv, dtype=str)

    sample_like_cols = [c for c in meta.columns if "sample" in c.lower() or "source name" in c.lower() or "title" in c.lower()]
    if not sample_like_cols:
        raise ValueError("No sample-like columns found in metadata CSV")

    meta = meta.copy()
    meta["matched_sample_number"] = pd.NA

    for col in sample_like_cols:
        nums = meta[col].map(sample_number_from_text)
        meta["matched_sample_number"] = meta["matched_sample_number"].fillna(nums)

    matched = meta.merge(
        reg[["sample_id", "sample_number"]],
        left_on="matched_sample_number",
        right_on="sample_number",
        how="left",
    )

    Path(args.outcsv).parent.mkdir(parents=True, exist_ok=True)
    matched.to_csv(args.outcsv, index=False)

    n_matched = matched["sample_id"].notna().sum()
    print("[ok] wrote matched metadata table:", args.outcsv)
    print("[info] rows:", len(matched))
    print("[info] matched rows:", int(n_matched))
    print("[info] unmatched rows:", int(len(matched) - n_matched))


if __name__ == "__main__":
    main()
