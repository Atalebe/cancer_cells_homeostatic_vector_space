from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


PLACEHOLDER_COLS = [
    "placeholder_condition",
    "placeholder_group",
    "placeholder_batch",
    "placeholder_note",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit GSE240704 sample annotations and detect placeholder-only columns."
    )
    parser.add_argument(
        "--sample-annotations",
        required=True,
        help="Path to sample annotations parquet or csv.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for annotation audit.",
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
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = read_table(args.sample_annotations)

    summary_rows = []
    for col in df.columns:
        s = df[col]
        non_null = int(s.notna().sum())
        unique_non_null = int(s.dropna().astype(str).nunique())
        top_values = s.dropna().astype(str).value_counts().head(10).to_dict()

        summary_rows.append(
            {
                "column": col,
                "dtype": str(s.dtype),
                "n_rows": int(len(df)),
                "non_null": non_null,
                "null": int(len(df) - non_null),
                "unique_non_null": unique_non_null,
                "all_missing": bool(non_null == 0),
                "single_value_non_null": bool(unique_non_null == 1 and non_null > 0),
                "top_values_json": json.dumps(top_values, ensure_ascii=False),
            }
        )

    audit = pd.DataFrame(summary_rows).sort_values(
        ["all_missing", "unique_non_null", "non_null"],
        ascending=[True, False, False],
    )
    audit.to_csv(outdir / "sample_annotations_audit.csv", index=False)

    placeholder_presence = []
    for col in PLACEHOLDER_COLS:
        if col in df.columns:
            vals = df[col].dropna().astype(str)
            placeholder_presence.append(
                {
                    "column": col,
                    "present": True,
                    "non_null": int(vals.shape[0]),
                    "unique_values": int(vals.nunique()),
                    "value_counts_json": json.dumps(vals.value_counts().head(20).to_dict(), ensure_ascii=False),
                }
            )
        else:
            placeholder_presence.append(
                {
                    "column": col,
                    "present": False,
                    "non_null": 0,
                    "unique_values": 0,
                    "value_counts_json": json.dumps({}, ensure_ascii=False),
                }
            )

    pd.DataFrame(placeholder_presence).to_csv(
        outdir / "placeholder_annotation_status.csv",
        index=False,
    )

    informative_cols = audit.loc[
        (~audit["all_missing"]) & (~audit["single_value_non_null"])
    ].copy()
    informative_cols.to_csv(outdir / "sample_annotations_informative_columns.csv", index=False)

    print("[ok] wrote annotation audit tables to", outdir)
    print("[info] total columns:", df.shape[1])
    print("[info] informative columns:", len(informative_cols))

    if "sample_id" in df.columns:
        print("[info] sample_id present")
    else:
        print("[warn] sample_id missing")


if __name__ == "__main__":
    main()
