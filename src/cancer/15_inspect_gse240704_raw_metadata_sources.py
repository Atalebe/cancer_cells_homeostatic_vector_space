from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


TEXT_EXTS = {".txt", ".tsv", ".csv", ".json", ".soft", ".idf", ".sdrf"}
MAX_PREVIEW_LINES = 80


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inspect raw metadata sources for GSE240704 and identify sample-level biological fields."
    )
    parser.add_argument(
        "--raw-dir",
        required=True,
        help="Raw data directory for GSE240704.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for metadata source inspection.",
    )
    return parser.parse_args()


def sniff_delimiter(header: str) -> str | None:
    if "\t" in header:
        return "\t"
    if "," in header:
        return ","
    if ";" in header:
        return ";"
    return None


def looks_sample_column(name: str) -> bool:
    x = name.strip().lower()
    return any(k in x for k in [
        "sample", "source name", "title", "characteristics", "patient",
        "group", "condition", "batch", "treatment", "cell", "status",
        "disease", "phenotype"
    ])


def text_preview(path: Path) -> list[str]:
    lines = []
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for i, line in enumerate(f):
                lines.append(line.rstrip("\n"))
                if i + 1 >= MAX_PREVIEW_LINES:
                    break
    except Exception as e:
        lines = [f"[error reading file: {e}]"]
    return lines


def inspect_tabular_file(path: Path) -> tuple[pd.DataFrame | None, dict]:
    preview = text_preview(path)
    nonempty = [x for x in preview if x.strip()]
    header = nonempty[0] if nonempty else ""
    delim = sniff_delimiter(header)

    meta = {
        "file": str(path),
        "delimiter_guess": delim,
        "n_preview_lines": len(preview),
        "header_preview": header[:500],
        "candidate_columns": [],
        "read_ok": False,
        "n_rows_preview": 0,
    }

    if delim is None:
        return None, meta

    try:
        df = pd.read_csv(path, sep=delim, nrows=200, dtype=str, low_memory=False)
        meta["read_ok"] = True
        meta["n_rows_preview"] = len(df)
        meta["candidate_columns"] = [c for c in df.columns if looks_sample_column(c)]
        return df, meta
    except Exception as e:
        meta["read_ok"] = False
        meta["read_error"] = str(e)
        return None, meta


def inspect_json_file(path: Path) -> dict:
    meta = {
        "file": str(path),
        "json_top_keys": [],
        "sample_like_keys": [],
        "read_ok": False,
    }
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            obj = json.load(f)
        if isinstance(obj, dict):
            keys = list(obj.keys())
            meta["json_top_keys"] = keys[:100]
            meta["sample_like_keys"] = [k for k in keys if looks_sample_column(k)]
        meta["read_ok"] = True
    except Exception as e:
        meta["read_ok"] = False
        meta["read_error"] = str(e)
    return meta


def main() -> None:
    args = parse_args()
    raw_dir = Path(args.raw_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    files = sorted([
        p for p in raw_dir.rglob("*")
        if p.is_file() and p.suffix.lower() in TEXT_EXTS
    ])

    source_rows = []
    candidate_tables = []

    for path in files:
        suffix = path.suffix.lower()

        if suffix == ".json":
            meta = inspect_json_file(path)
            source_rows.append(meta)
            continue

        df, meta = inspect_tabular_file(path)
        source_rows.append(meta)

        if df is not None and meta["candidate_columns"]:
            keep = [c for c in df.columns if c in meta["candidate_columns"]]
            cand = df[keep].copy()
            cand.insert(0, "source_file", str(path))
            candidate_tables.append(cand)

    pd.DataFrame(source_rows).to_csv(
        outdir / "raw_metadata_source_inventory.csv",
        index=False,
    )

    if candidate_tables:
        union = pd.concat(candidate_tables, ignore_index=True, sort=False)
        union.to_csv(outdir / "raw_metadata_candidate_columns_union.csv", index=False)
        print("[ok] wrote candidate metadata union:", outdir / "raw_metadata_candidate_columns_union.csv")
    else:
        print("[warn] no candidate tabular metadata fields found")

    print("[ok] wrote metadata source inventory:", outdir / "raw_metadata_source_inventory.csv")
    print("[info] files inspected:", len(files))


if __name__ == "__main__":
    main()
