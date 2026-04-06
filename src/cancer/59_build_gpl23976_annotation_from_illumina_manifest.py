#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


DEFAULT_USECOLS = [
    "IlmnID",
    "Name",
    "Infinium_Design_Type",
    "Color_Channel",
    "Genome_Build",
    "CHR",
    "MAPINFO",
    "Strand",
    "UCSC_RefGene_Name",
    "UCSC_RefGene_Accession",
    "UCSC_RefGene_Group",
    "UCSC_CpG_Islands_Name",
    "Relation_to_UCSC_CpG_Island",
    "Phantom4_Enhancers",
    "HMM_Island",
    "Regulatory_Feature_Name",
    "Regulatory_Feature_Group",
    "DNase_Hypersensitivity_NAME",
    "OpenChromatin_NAME",
    "TFBS_NAME",
]


def detect_assay_header_line(path: Path, max_lines: int = 200) -> int:
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for i in range(max_lines):
            line = f.readline()
            if not line:
                break
            if line.startswith("IlmnID,Name,"):
                return i
    raise RuntimeError("could not find assay header line beginning with 'IlmnID,Name,'")


def read_manifest_csv(path: Path, usecols: list[str]) -> pd.DataFrame:
    header_line = detect_assay_header_line(path)
    df = pd.read_csv(
        path,
        skiprows=header_line,
        usecols=lambda c: c in usecols,
        low_memory=False,
    )
    return df


def normalize_manifest(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    rename_map = {
        "IlmnID": "probe_id",
        "Name": "ID_REF",
        "Infinium_Design_Type": "design_type",
        "Color_Channel": "color_channel",
        "Genome_Build": "genome_build",
        "CHR": "chromosome",
        "MAPINFO": "mapinfo",
        "Strand": "strand",
        "UCSC_RefGene_Name": "ucsc_refgene_name",
        "UCSC_RefGene_Accession": "ucsc_refgene_accession",
        "UCSC_RefGene_Group": "refgene_group",
        "UCSC_CpG_Islands_Name": "cpg_island_name",
        "Relation_to_UCSC_CpG_Island": "relation_to_cpg_island",
        "Phantom4_Enhancers": "phantom4_enhancers",
        "HMM_Island": "hmm_island",
        "Regulatory_Feature_Name": "regulatory_feature_name",
        "Regulatory_Feature_Group": "regulatory_feature_group",
        "DNase_Hypersensitivity_NAME": "dnase_name",
        "OpenChromatin_NAME": "openchromatin_name",
        "TFBS_NAME": "tfbs_name",
    }
    out = out.rename(columns=rename_map)

    for col in out.columns:
        if out[col].dtype == object:
            out[col] = out[col].astype("string").str.strip()

    if "ID_REF" not in out.columns and "probe_id" in out.columns:
        out["ID_REF"] = out["probe_id"]

    # Keep original Illumina names too, if present, because they help later debugging
    inverse_keep = {
        "ucsc_refgene_name": "UCSC_RefGene_Name",
        "ucsc_refgene_accession": "UCSC_RefGene_Accession",
        "refgene_group": "UCSC_RefGene_Group",
        "cpg_island_name": "UCSC_CpG_Islands_Name",
        "relation_to_cpg_island": "Relation_to_UCSC_CpG_Island",
    }
    for norm_col, raw_col in inverse_keep.items():
        if norm_col in out.columns and raw_col not in out.columns:
            out[raw_col] = out[norm_col]

    out = out.dropna(subset=["ID_REF"])
    out = out[out["ID_REF"].astype("string").str.len() > 0]
    out = out.drop_duplicates(subset=["ID_REF"], keep="first").reset_index(drop=True)

    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-file", required=True)
    ap.add_argument("--outcsv", required=True)
    ap.add_argument("--outparquet", required=True)
    ap.add_argument("--summary-json", required=True)
    ap.add_argument("--head-csv", required=True)
    args = ap.parse_args()

    input_path = Path(args.input_file)
    outcsv = Path(args.outcsv)
    outparquet = Path(args.outparquet)
    summary_json = Path(args.summary_json)
    head_csv = Path(args.head_csv)

    outcsv.parent.mkdir(parents=True, exist_ok=True)

    raw_df = read_manifest_csv(input_path, DEFAULT_USECOLS)
    norm_df = normalize_manifest(raw_df)

    norm_df.to_csv(outcsv, index=False)
    norm_df.to_parquet(outparquet, index=False)
    norm_df.head(50).to_csv(head_csv, index=False)

    summary = {
        "input_file": str(input_path),
        "n_rows_loaded": int(len(raw_df)),
        "n_rows_normalized": int(len(norm_df)),
        "n_unique_id_ref": int(norm_df["ID_REF"].nunique()) if "ID_REF" in norm_df.columns else 0,
        "columns_written": list(norm_df.columns),
        "usecols": DEFAULT_USECOLS,
    }
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[ok] wrote normalized manifest csv: {outcsv}")
    print(f"[ok] wrote normalized manifest parquet: {outparquet}")
    print(f"[ok] wrote summary: {summary_json}")
    print("[info] summary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
