#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


ANNOTATION_HINTS = [
    "gene", "symbol", "chr", "chrom", "feature", "island",
    "strand", "mapinfo", "position", "accession", "refgene",
    "relation", "region", "name", "description", "transcript",
    "ucsc", "ilmn", "probe", "id"
]

PROBE_HINTS = ["probe", "id", "ilmn", "cg", "spot"]
GENE_HINTS = ["gene", "symbol", "refgene", "name"]
FEATURE_HINTS = ["feature", "relation", "group", "class", "region"]
ISLAND_HINTS = ["island"]
CHR_HINTS = ["chr", "chrom", "mapinfo", "position"]


def safe_read_text(path: Path, n_bytes: int = 50000) -> str:
    with open(path, "rb") as f:
        raw = f.read(n_bytes)
    return raw.decode("utf-8", errors="replace")


def looks_binary(text: str) -> bool:
    return "\x00" in text


def detect_delim_from_text(text: str) -> Optional[str]:
    lines = [ln for ln in text.splitlines()[:20] if ln.strip()]
    if not lines:
        return None
    tab_votes = sum("\t" in ln for ln in lines)
    comma_votes = sum("," in ln for ln in lines)
    if tab_votes >= 2:
        return "\t"
    if comma_votes >= 2:
        return ","
    return None


def score_columns(cols: List[str]) -> int:
    score = 0
    for c in cols:
        cl = c.lower()
        for hint in ANNOTATION_HINTS:
            if hint in cl:
                score += 1
    return score


def pick_candidate_columns(cols: List[str], hints: List[str]) -> List[str]:
    out = []
    for c in cols:
        cl = c.lower()
        if any(h in cl for h in hints):
            out.append(c)
    return out


def inspect_df(df: pd.DataFrame, source_name: str) -> Dict[str, object]:
    cols = [str(c) for c in df.columns]
    return {
        "source_name": source_name,
        "n_rows_preview": int(df.shape[0]),
        "n_cols_preview": int(df.shape[1]),
        "columns_preview": cols[:100],
        "annotation_col_score": int(score_columns(cols)),
        "candidate_probe_cols": pick_candidate_columns(cols, PROBE_HINTS),
        "candidate_gene_cols": pick_candidate_columns(cols, GENE_HINTS),
        "candidate_feature_cols": pick_candidate_columns(cols, FEATURE_HINTS),
        "candidate_island_cols": pick_candidate_columns(cols, ISLAND_HINTS),
        "candidate_chr_cols": pick_candidate_columns(cols, CHR_HINTS),
    }


def preview_soft_or_text(path: Path) -> Dict[str, object]:
    header_lines = []
    table_lines = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f):
            s = line.rstrip("\n")
            if s.startswith("!platform_table_begin") or s.startswith("!Platform_table_begin"):
                table_lines.append(s)
            elif s.startswith("!platform_table_end") or s.startswith("!Platform_table_end"):
                table_lines.append(s)
            elif s.startswith("#") or s.startswith("!"):
                header_lines.append(s)
            elif "\t" in s and len(table_lines) < 20:
                table_lines.append(s)
            if i > 400:
                break

    return {
        "mode": "soft_or_text_preview",
        "read_ok": False,
        "source_name": None,
        "n_rows_preview": None,
        "n_cols_preview": None,
        "annotation_col_score": None,
        "candidate_probe_cols": [],
        "candidate_gene_cols": [],
        "candidate_feature_cols": [],
        "candidate_island_cols": [],
        "candidate_chr_cols": [],
        "columns_preview": [],
        "notes": json.dumps({
            "header_preview": header_lines[:20],
            "table_preview": table_lines[:20],
        }),
    }


def try_html_tables(path: Path) -> Optional[Dict[str, object]]:
    try:
        tables = pd.read_html(str(path))
    except Exception:
        return None

    if not tables:
        return None

    best_info = None
    best_score = -1
    best_idx = -1
    for i, df in enumerate(tables):
        info = inspect_df(df.astype(str), f"html_table_{i}")
        if info["annotation_col_score"] > best_score:
            best_score = info["annotation_col_score"]
            best_info = info
            best_idx = i

    if best_info is None:
        return None

    best_info["mode"] = "html_table"
    best_info["read_ok"] = True
    best_info["source_name"] = f"html_table_{best_idx}"
    best_info["notes"] = ""
    return best_info


def try_delim_table(path: Path, text: str) -> Optional[Dict[str, object]]:
    delim = detect_delim_from_text(text)
    if delim is None:
        return None
    try:
        df = pd.read_csv(path, sep=delim, nrows=200, dtype=str, low_memory=False)
    except Exception:
        return None
    info = inspect_df(df, "detected_delim_table")
    info["mode"] = "detected_delim_table"
    info["read_ok"] = True
    info["notes"] = ""
    return info


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows: List[Dict[str, object]] = []
    union_cols: List[Dict[str, str]] = []
    best_candidate: Optional[Dict[str, object]] = None
    best_score = -1

    for path in sorted(indir.rglob("*")):
        if not path.is_file():
            continue

        row: Dict[str, object] = {
            "file": str(path),
            "suffix": path.suffix.lower(),
            "size_bytes": path.stat().st_size,
            "mode": None,
            "read_ok": False,
            "source_name": None,
            "n_rows_preview": None,
            "n_cols_preview": None,
            "annotation_col_score": None,
            "candidate_probe_cols": [],
            "candidate_gene_cols": [],
            "candidate_feature_cols": [],
            "candidate_island_cols": [],
            "candidate_chr_cols": [],
            "columns_preview": [],
            "notes": "",
        }

        try:
            text = safe_read_text(path)

            if looks_binary(text):
                row["notes"] = "binary_or_non_text"
                rows.append(row)
                continue

            lower = text.lower()

            # Important: HTML first
            if "<html" in lower or "<table" in lower or path.suffix.lower() in {".html", ".htm", ".cgi"}:
                info = try_html_tables(path)
                if info is None:
                    row["mode"] = "html_unparsed"
                    row["notes"] = "html_detected_but_no_table_parsed"
                else:
                    row.update(info)

            else:
                info = try_delim_table(path, text)
                if info is not None:
                    row.update(info)
                else:
                    row.update(preview_soft_or_text(path))

            if row["read_ok"]:
                for c in row["columns_preview"]:
                    union_cols.append({"file": str(path), "column": c})
                score = int(row["annotation_col_score"])
                if score > best_score:
                    best_score = score
                    best_candidate = {
                        "file": str(path),
                        "mode": row["mode"],
                        "source_name": row["source_name"],
                        "score": score,
                        "columns_preview": row["columns_preview"],
                    }

        except Exception as e:
            row["notes"] = repr(e)

        rows.append(row)

    pd.DataFrame(rows).to_csv(outdir / "bundle_file_inventory.csv", index=False)
    pd.DataFrame(union_cols).to_csv(outdir / "bundle_column_union.csv", index=False)

    summary = {
        "indir": str(indir),
        "n_files": len(rows),
        "best_candidate_file": None if best_candidate is None else best_candidate["file"],
        "best_candidate_mode": None if best_candidate is None else best_candidate["mode"],
        "best_candidate_source_name": None if best_candidate is None else best_candidate["source_name"],
        "best_candidate_score": best_score,
        "best_candidate_columns_preview": None if best_candidate is None else best_candidate["columns_preview"],
    }

    (outdir / "bundle_inspection_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print(f"[ok] wrote bundle inspection to {outdir}")
    print(f"[info] files inspected: {len(rows)}")
    print(f"[info] best candidate: {summary['best_candidate_file']}")
    print(f"[info] best mode: {summary['best_candidate_mode']}")
    print(f"[info] best score: {best_score}")


if __name__ == "__main__":
    main()
