#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd


CANDIDATE_ANNOTATION_FILES = [
    "data/raw/gse240704/annotation.csv",
    "data/raw/gse240704/annotation.tsv",
    "data/raw/gse240704/annotation.txt",
    "data/raw/gse240704/IlluminaHumanMethylationEPIC_annotation.csv",
    "data/raw/gse240704/IlluminaHumanMethylation450k_annotation.csv",
    "data/reference/illumina_methylation_annotation.csv",
    "data/reference/illumina_methylation_annotation.tsv",
]


def first_existing(paths: List[str]) -> Path | None:
    for p in paths:
        q = Path(p)
        if q.exists():
            return q
    return None


def read_table_auto(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    return pd.read_csv(path, sep=None, engine="python")


def normalize_cols(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [str(c).strip() for c in out.columns]
    return out


def choose_probe_col(df: pd.DataFrame) -> str:
    preferred = [
        "ID_REF",
        "IlmnID",
        "Probe_ID",
        "probe_id",
        "Name",
        "TargetID",
    ]
    for c in preferred:
        if c in df.columns:
            return c
    for c in df.columns:
        low = c.lower()
        if "probe" in low or "ilmn" in low or "id_ref" in low or low == "name":
            return c
    raise ValueError("Could not identify probe ID column in annotation file")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contrast-csv", required=True)
    ap.add_argument("--outcsv", required=True)
    ap.add_argument("--outjson", required=True)
    ap.add_argument("--annotation-file", default=None)
    args = ap.parse_args()

    contrast = normalize_cols(pd.read_csv(args.contrast_csv))

    ann_path = Path(args.annotation_file) if args.annotation_file else first_existing(CANDIDATE_ANNOTATION_FILES)

    summary = {
        "contrast_csv": args.contrast_csv,
        "annotation_file_used": None,
        "annotation_rows": 0,
        "annotation_probe_col": None,
        "n_contrast_rows": int(len(contrast)),
        "n_annotated_rows": 0,
    }

    if ann_path is None:
        contrast["direction"] = np.where(
            contrast["delta_median_a_minus_b"] >= 0, "higher_in_a", "lower_in_a"
        )
        contrast.to_csv(args.outcsv, index=False)
        with open(args.outjson, "w") as f:
            json.dump(summary, f, indent=2)
        print("[warn] no annotation file found, wrote unannotated contrast table")
        print(f"[ok] wrote: {args.outcsv}")
        print(f"[ok] wrote: {args.outjson}")
        return

    ann = normalize_cols(read_table_auto(ann_path))
    probe_col = choose_probe_col(ann)

    summary["annotation_file_used"] = str(ann_path)
    summary["annotation_rows"] = int(len(ann))
    summary["annotation_probe_col"] = probe_col

    ann = ann.rename(columns={probe_col: "ID_REF"})
    keep_cols = ["ID_REF"] + [
        c for c in ann.columns
        if c != "ID_REF" and (
            "gene" in c.lower()
            or "feature" in c.lower()
            or "island" in c.lower()
            or "chr" in c.lower()
            or "map" in c.lower()
            or "relation" in c.lower()
            or "regulatory" in c.lower()
        )
    ]
    ann_keep = ann[keep_cols].drop_duplicates(subset=["ID_REF"])

    merged = contrast.merge(ann_keep, on="ID_REF", how="left")
    merged["direction"] = np.where(
        merged["delta_median_a_minus_b"] >= 0, "higher_in_a", "lower_in_a"
    )

    annotation_cols = [c for c in merged.columns if c not in contrast.columns and c != "direction"]
    if annotation_cols:
        summary["n_annotated_rows"] = int(merged[annotation_cols].notna().any(axis=1).sum())

    Path(args.outcsv).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.outcsv, index=False)
    with open(args.outjson, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[ok] wrote annotated contrast: {args.outcsv}")
    print(f"[ok] wrote summary: {args.outjson}")
    print(f"[info] annotation file: {ann_path}")
    print(f"[info] annotated rows: {summary['n_annotated_rows']}")
    print(f"[info] annotation columns kept: {len(annotation_cols)}")


if __name__ == "__main__":
    main()
