from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract sample-level metadata from a GEO series matrix text file."
    )
    parser.add_argument("--series-matrix", required=True)
    parser.add_argument("--outcsv", required=True)
    return parser.parse_args()


def clean_field(x: str) -> str:
    x = str(x).strip()
    if x.startswith('"') and x.endswith('"'):
        x = x[1:-1]
    return x.strip()


def main() -> None:
    args = parse_args()
    matrix_path = Path(args.series_matrix)

    rows = []
    with open(matrix_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("!Sample_"):
                continue
            line = line.rstrip("\n")
            parsed = next(csv.reader([line], delimiter="\t", quotechar='"'))
            key = parsed[0].strip()
            values = [clean_field(x) for x in parsed[1:]]
            rows.append((key, values))

    if not rows:
        raise ValueError("No !Sample_ lines found in series matrix.")

    sample_count = max(len(v) for _, v in rows)
    sample_ids = [f"SAMPLE {i}" for i in range(1, sample_count + 1)]

    data = {"sample_id": sample_ids}
    for key, values in rows:
        col = key.replace("!Sample_", "").strip()
        padded = values + [""] * (sample_count - len(values))
        data[col] = padded[:sample_count]

    df = pd.DataFrame(data)
    Path(args.outcsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.outcsv, index=False)

    print("[ok] wrote extracted series-matrix metadata:", args.outcsv)
    print("[info] rows:", len(df))
    print("[info] columns:", len(df.columns))
    print("[info] first columns:")
    for c in df.columns[:20]:
        print(" ", c)


if __name__ == "__main__":
    main()
