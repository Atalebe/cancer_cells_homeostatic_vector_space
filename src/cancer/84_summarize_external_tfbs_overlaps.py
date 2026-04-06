#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
from collections import Counter

import pandas as pd


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
BASE = REPO_ROOT / "results" / "gse240704" / "external_tfbs_overlap"
REMAP_DIR = BASE / "remap_json"
UCSC_DIR = BASE / "ucsc_json"
OUTDIR = BASE / "summary"
OUTDIR.mkdir(parents=True, exist_ok=True)


def load_json(path: Path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def flatten_records(obj):
    if isinstance(obj, list):
        return obj
    if isinstance(obj, dict):
        # common wrappers
        for key in ["items", "results", "data", "peaks", "track", "wgEncodeRegTfbsClusteredV3"]:
            if key in obj and isinstance(obj[key], list):
                return obj[key]
        return [obj]
    return []


def guess_tf_fields(df: pd.DataFrame) -> list[str]:
    candidates = []
    for c in df.columns:
        lc = c.lower()
        if any(k in lc for k in ["tf", "factor", "antigen", "target", "name"]):
            candidates.append(c)
    return candidates


def summarize_remap_file(path: Path) -> dict:
    try:
        obj = load_json(path)
    except Exception as e:
        return {"label": path.stem, "source": "remap", "status": f"failed_to_load_json: {e}"}

    records = flatten_records(obj)
    if not records:
        return {"label": path.stem, "source": "remap", "status": "no_records"}

    df = pd.json_normalize(records)
    tf_fields = guess_tf_fields(df)

    tf_counter = Counter()
    if tf_fields:
        for col in tf_fields:
            for val in df[col].dropna().astype(str):
                if val.strip():
                    tf_counter[val.strip()] += 1

    return {
        "label": path.stem,
        "source": "remap",
        "status": "ok",
        "n_rows": int(len(df)),
        "candidate_tf_fields": tf_fields,
        "top_tf_like_values": dict(tf_counter.most_common(20)),
        "all_columns": list(df.columns),
    }


def summarize_ucsc_file(path: Path) -> dict:
    try:
        obj = load_json(path)
    except Exception as e:
        return {"label": path.stem, "source": "ucsc", "status": f"failed_to_load_json: {e}"}

    records = flatten_records(obj)
    if not records:
        return {"label": path.stem, "source": "ucsc", "status": "no_records"}

    df = pd.json_normalize(records)
    tf_fields = guess_tf_fields(df)

    tf_counter = Counter()
    if tf_fields:
        for col in tf_fields:
            for val in df[col].dropna().astype(str):
                if val.strip():
                    tf_counter[val.strip()] += 1

    return {
        "label": path.stem,
        "source": "ucsc",
        "status": "ok",
        "n_rows": int(len(df)),
        "candidate_tf_fields": tf_fields,
        "top_tf_like_values": dict(tf_counter.most_common(20)),
        "all_columns": list(df.columns),
    }


def main():
    rows = []

    for path in sorted(REMAP_DIR.glob("*.json")):
        rows.append(summarize_remap_file(path))

    for path in sorted(UCSC_DIR.glob("*.json")):
        rows.append(summarize_ucsc_file(path))

    summary_df = pd.DataFrame(rows)
    summary_csv = OUTDIR / "external_tfbs_overlap_summary.csv"
    summary_df.to_csv(summary_csv, index=False)

    summary_json = OUTDIR / "external_tfbs_overlap_summary.json"
    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=2)

    print(summary_df)
    print(f"[ok] wrote {summary_csv}")
    print(f"[ok] wrote {summary_json}")


if __name__ == "__main__":
    main()
