from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build editable metadata mapping template for GSE240704 samples."
    )
    parser.add_argument(
        "--sample-registry",
        required=True,
        help="Path to sample registry csv.",
    )
    parser.add_argument(
        "--sample-annotations",
        required=False,
        default=None,
        help="Optional existing sample annotations parquet or csv.",
    )
    parser.add_argument(
        "--outcsv",
        required=True,
        help="Output CSV path for manual metadata mapping.",
    )
    return parser.parse_args()


def read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() == ".parquet":
        return pd.read_parquet(p)
    if p.suffix.lower() == ".csv":
        return pd.read_csv(p)
    raise ValueError(f"Unsupported file format: {p}")


def main() -> None:
    args = parse_args()

    reg = read_table(args.sample_registry)
    if "sample_id" not in reg.columns:
        raise ValueError("sample_registry must contain sample_id")

    cols = [c for c in ["sample_id", "sample_number"] if c in reg.columns]
    out = reg[cols].copy().sort_values(cols[0])

    if args.sample_annotations:
        ann = read_table(args.sample_annotations)
        keep = [c for c in ann.columns if c in ["sample_id", "placeholder_condition", "placeholder_group", "placeholder_batch", "placeholder_note"]]
        if "sample_id" in keep:
            out = out.merge(ann[keep], on="sample_id", how="left")

    for col in [
        "biological_condition",
        "disease_group",
        "cell_type",
        "tumor_status",
        "batch",
        "patient_id",
        "timepoint",
        "sex",
        "age",
        "source_note",
    ]:
        if col not in out.columns:
            out[col] = pd.NA

    Path(args.outcsv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.outcsv, index=False)

    print("[ok] wrote editable metadata mapping template:", args.outcsv)
    print("[info] rows:", len(out))


if __name__ == "__main__":
    main()
