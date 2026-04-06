#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--soft-gz", required=True)
    ap.add_argument("--outcsv", required=True)
    ap.add_argument("--outparquet", required=True)
    ap.add_argument("--outjson", required=True)
    args = ap.parse_args()

    soft_path = Path(args.soft_gz)
    if not soft_path.exists():
        raise SystemExit(f"missing input: {soft_path}")

    in_platform = False
    saw_header = False
    header = None
    rows = []

    with gzip.open(soft_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")

            if line.startswith("^PLATFORM"):
                in_platform = True
                saw_header = False
                header = None
                continue

            if in_platform and line.startswith("^") and not line.startswith("^PLATFORM"):
                break

            if not in_platform:
                continue

            if line.startswith("!platform_table_begin"):
                saw_header = False
                header = None
                continue

            if line.startswith("!platform_table_end"):
                break

            if line.startswith("#"):
                continue

            if line.startswith("!"):
                continue

            if not saw_header:
                header = line.split("\t")
                saw_header = True
                continue

            vals = line.split("\t")
            if header is None:
                continue

            if len(vals) < len(header):
                vals = vals + [""] * (len(header) - len(vals))
            elif len(vals) > len(header):
                vals = vals[: len(header)]

            rows.append(vals)

    if header is None or len(rows) == 0:
        raise SystemExit("failed to parse platform annotation table from SOFT file")

    df = pd.DataFrame(rows, columns=header)

    outcsv = Path(args.outcsv)
    outparquet = Path(args.outparquet)
    outjson = Path(args.outjson)
    outcsv.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(outcsv, index=False)
    df.to_parquet(outparquet, index=False)

    summary = {
        "soft_gz": str(soft_path),
        "n_rows": int(len(df)),
        "n_columns": int(len(df.columns)),
        "columns": list(map(str, df.columns)),
    }

    with open(outjson, "w") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote platform annotation csv:", outcsv)
    print("[ok] wrote platform annotation parquet:", outparquet)
    print("[ok] wrote summary:", outjson)
    print("[info] rows:", len(df))
    print("[info] columns:", len(df.columns))


if __name__ == "__main__":
    main()
