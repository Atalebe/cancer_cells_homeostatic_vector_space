from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


PREFERRED_COL_HINTS = [
    "sample",
    "source name",
    "title",
    "characteristics",
    "patient",
    "group",
    "condition",
    "batch",
    "treatment",
    "cell",
    "status",
    "disease",
    "phenotype",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract a sample metadata candidate table from a chosen raw metadata file."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to selected metadata source file, usually SDRF/TSV/CSV.",
    )
    parser.add_argument(
        "--outcsv",
        required=True,
        help="Output CSV path.",
    )
    return parser.parse_args()


def sniff_delimiter(path: Path) -> str:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        first = f.readline()
    if "\t" in first:
        return "\t"
    if "," in first:
        return ","
    if ";" in first:
        return ";"
    raise ValueError(f"Could not determine delimiter for {path}")


def keep_col(name: str) -> bool:
    x = name.strip().lower()
    return any(k in x for k in PREFERRED_COL_HINTS)


def main() -> None:
    args = parse_args()
    path = Path(args.input)
    sep = sniff_delimiter(path)

    df = pd.read_csv(path, sep=sep, dtype=str, low_memory=False)
    cols = [c for c in df.columns if keep_col(c)]

    if not cols:
        raise ValueError("No likely sample metadata columns found in selected file.")

    out = df[cols].copy()
    out.insert(0, "source_file", str(path))

    Path(args.outcsv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.outcsv, index=False)

    print("[ok] wrote extracted sample metadata candidate table:", args.outcsv)
    print("[info] rows:", len(out))
    print("[info] columns:")
    for c in out.columns:
        print(" ", c)


if __name__ == "__main__":
    main()
