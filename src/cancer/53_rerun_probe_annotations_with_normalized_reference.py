#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ANNOTATION_PRIORITY = [
    "ID_REF",
    "probe_id",
    "ucsc_refgene_name",
    "ucsc_refgene_accession",
    "refgene_group",
    "cpg_island_name",
    "relation_to_cpg_island",
    "regulatory_feature_group",
    "regulatory_feature_name",
    "phantom4_enhancers",
    "hmm_island",
    "chromosome",
    "mapinfo",
    "strand",
    "design_type",
    "color_channel",
    "UCSC_RefGene_Name",
    "UCSC_RefGene_Accession",
    "UCSC_RefGene_Group",
    "UCSC_CpG_Islands_Name",
    "Relation_to_UCSC_CpG_Island",
]


def choose_probe_col(df: pd.DataFrame) -> str | None:
    for c in ["ID_REF", "probe_id", "Name", "ID", "IlmnID"]:
        if c in df.columns:
            return c
    return None


def reannotate_one(annotation_df: pd.DataFrame, input_csv: str, outdir: Path) -> dict:
    df = pd.read_csv(input_csv)
    probe_col = choose_probe_col(df)
    if probe_col is None:
        raise RuntimeError(f"could not determine probe column in {input_csv}")

    ann_probe_col = choose_probe_col(annotation_df)
    if ann_probe_col is None:
        raise RuntimeError("could not determine probe column in annotation parquet")

    ann_cols = [c for c in ANNOTATION_PRIORITY if c in annotation_df.columns]
    ann = annotation_df[ann_cols].copy()

    merged = df.merge(
        ann,
        left_on=probe_col,
        right_on=ann_probe_col,
        how="left",
        suffixes=("", "_ann"),
    )

    stem = Path(input_csv).stem
    out_csv = outdir / f"{stem}_normalized_annotation.csv"
    out_json = outdir / f"{stem}_normalized_annotation_summary.json"

    merged.to_csv(out_csv, index=False)

    summary = {
        "input_csv": input_csv,
        "probe_col_input": probe_col,
        "probe_col_annotation": ann_probe_col,
        "n_input_rows": int(len(df)),
        "n_output_rows": int(len(merged)),
        "annotation_columns_used": ann_cols,
        "n_annotated_rows": int(merged[ann_probe_col].notna().sum()) if ann_probe_col in merged.columns else 0,
        "output_csv": str(out_csv),
    }
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[ok] wrote: {out_csv}")
    print(f"[ok] wrote: {out_json}")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotation-parquet", required=True)
    ap.add_argument("--contrast-csvs", nargs="+", required=True)
    ap.add_argument("--s-axis-driver-csv", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    annotation_df = pd.read_parquet(args.annotation_parquet)

    batch = []
    for p in args.contrast_csvs + [args.s_axis_driver_csv]:
        batch.append(reannotate_one(annotation_df, p, outdir))

    batch_json = outdir / "reannotation_batch_summary.json"
    with open(batch_json, "w") as f:
        json.dump(batch, f, indent=2)

    print(f"[ok] wrote batch summary: {batch_json}")


if __name__ == "__main__":
    main()
