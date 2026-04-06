#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def choose_probe_col(df: pd.DataFrame) -> str | None:
    preferred = [
        "ID_REF",
        "probe_id",
        "Probe_ID",
        "IlmnID",
        "ilmnid",
        "Name",
        "name",
        "cg_id",
        "CpG",
        "TargetID",
    ]
    for col in preferred:
        if col in df.columns:
            return col

    for col in df.columns:
        s = df[col].dropna().astype(str)
        if s.empty:
            continue
        frac = s.head(500).str.match(r"^cg\d+$", na=False).mean()
        if frac > 0.20:
            return col
    return None


def choose_annotation_columns(df: pd.DataFrame):
    wanted = []
    priority = [
        "ID_REF",
        "Name",
        "gene",
        "Gene",
        "Gene_Symbol",
        "gene_symbol",
        "UCSC_RefGene_Name",
        "RefGene_Name",
        "UCSC_RefGene_Group",
        "Relation_to_Island",
        "Island",
        "CHR",
        "MAPINFO",
        "Chromosome",
        "Coordinate",
        "Feature",
        "Description",
    ]
    for col in priority:
        if col in df.columns and col not in wanted:
            wanted.append(col)

    for col in df.columns:
        if col not in wanted and df[col].dtype == "object":
            wanted.append(col)
        if len(wanted) >= 20:
            break
    return wanted


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contrast-csv", required=True)
    ap.add_argument("--annotation-parquet", required=True)
    ap.add_argument("--outcsv", required=True)
    ap.add_argument("--outjson", required=True)
    ap.add_argument("--probe-col-contrast", default="ID_REF")
    ap.add_argument("--probe-col-annotation", default=None)
    args = ap.parse_args()

    contrast = pd.read_csv(args.contrast_csv)
    ann = pd.read_parquet(args.annotation_parquet)

    probe_col_contrast = args.probe_col_contrast
    if probe_col_contrast not in contrast.columns:
        raise SystemExit(f"contrast probe column missing: {probe_col_contrast}")

    probe_col_annotation = args.probe_col_annotation or choose_probe_col(ann)
    if probe_col_annotation is None:
        raise SystemExit("could not determine probe column in annotation parquet")

    ann_keep = choose_annotation_columns(ann)
    if probe_col_annotation not in ann_keep:
        ann_keep = [probe_col_annotation] + ann_keep

    ann2 = ann[ann_keep].copy()
    ann2[probe_col_annotation] = ann2[probe_col_annotation].astype(str)
    contrast[probe_col_contrast] = contrast[probe_col_contrast].astype(str)

    ann2 = ann2.drop_duplicates(subset=[probe_col_annotation])

    merged = contrast.merge(
        ann2,
        left_on=probe_col_contrast,
        right_on=probe_col_annotation,
        how="left",
    )

    outcsv = Path(args.outcsv)
    outcsv.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(outcsv, index=False)

    n_annotated = int(merged[probe_col_annotation].notna().sum())

    summary = {
        "contrast_csv": args.contrast_csv,
        "annotation_parquet": args.annotation_parquet,
        "probe_col_contrast": probe_col_contrast,
        "probe_col_annotation": probe_col_annotation,
        "n_contrast_rows": int(len(contrast)),
        "n_annotation_rows": int(len(ann)),
        "n_annotated_rows": n_annotated,
        "annotation_columns_used": ann_keep,
    }

    outjson = Path(args.outjson)
    with open(outjson, "w") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote annotated contrast:", outcsv)
    print("[ok] wrote summary:", outjson)
    print("[info] probe_col_annotation:", probe_col_annotation)
    print("[info] n_annotated_rows:", n_annotated)


if __name__ == "__main__":
    main()
