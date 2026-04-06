#!/usr/bin/env python3

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Optional

import pandas as pd
import yaml


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"

CANDIDATE_ANNOTATION_FILES = [
    REPO_ROOT / "data" / "reference" / "gencode" / "gencode_v38_gene_annotation.tsv",
    REPO_ROOT / "data" / "reference" / "gencode" / "gencode_gene_annotation.tsv",
    REPO_ROOT / "data" / "reference" / "annotations" / "ensembl_gene_annotation.tsv",
    REPO_ROOT / "data" / "reference" / "annotations" / "gene_annotation.tsv",
]


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def strip_version(gene_id: str) -> str:
    return str(gene_id).split(".", 1)[0]


def find_annotation_file() -> Optional[Path]:
    for path in CANDIDATE_ANNOTATION_FILES:
        if path.exists():
            return path
    return None


def classify_symbol(symbol: str) -> dict:
    s = "" if pd.isna(symbol) else str(symbol)
    return {
        "is_mito": bool(re.match(r"^MT-", s)),
        "is_ribo": bool(re.match(r"^RPL|^RPS|^MRPL|^MRPS", s)),
        "is_hla": bool(re.match(r"^HLA-", s)),
    }


def annotate_table(df: pd.DataFrame, ann: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["gene_id_stripped"] = out["gene_id"].astype(str).map(strip_version)
    ann["gene_id_stripped"] = ann["gene_id"].astype(str).map(strip_version)

    merged = out.merge(
        ann[["gene_id_stripped", "gene_symbol", "gene_biotype"]],
        on="gene_id_stripped",
        how="left",
    )

    flags = merged["gene_symbol"].map(classify_symbol).apply(pd.Series)
    merged = pd.concat([merged, flags], axis=1)
    return merged


def summarize_table(df: pd.DataFrame, label: str) -> dict:
    return {
        "label": label,
        "n_rows": int(len(df)),
        "n_with_symbol": int(df["gene_symbol"].notna().sum()),
        "n_mito": int(df["is_mito"].fillna(False).sum()),
        "n_ribo": int(df["is_ribo"].fillna(False).sum()),
        "n_hla": int(df["is_hla"].fillna(False).sum()),
        "top_gene_symbols_first20": df["gene_symbol"].fillna(df["gene_id"]).head(20).tolist(),
    }


def main() -> None:
    cfg = read_config()
    in_dir = REPO_ROOT / cfg["results_dir"] / "d1_d2_marker_scan_clean"
    out_dir = REPO_ROOT / cfg["results_dir"] / "d1_d2_marker_annotation"
    out_dir.mkdir(parents=True, exist_ok=True)

    ann_path = find_annotation_file()
    if ann_path is None:
        raise FileNotFoundError(
            "No annotation TSV found. Put a file with columns gene_id, gene_symbol, gene_biotype in one of: "
            + ", ".join(str(p) for p in CANDIDATE_ANNOTATION_FILES)
        )

    ann = pd.read_csv(ann_path, sep="\t")
    required = {"gene_id", "gene_symbol", "gene_biotype"}
    missing = required - set(ann.columns)
    if missing:
        raise ValueError(f"Annotation file missing columns: {sorted(missing)}")

    d1_high = pd.read_csv(in_dir / "d1_high_top200_clean.csv")
    d1_low = pd.read_csv(in_dir / "d1_low_top200_clean.csv")

    d1_high_ann = annotate_table(d1_high, ann)
    d1_low_ann = annotate_table(d1_low, ann)

    d1_high_ann.to_csv(out_dir / "d1_high_top200_clean_annotated.csv", index=False)
    d1_low_ann.to_csv(out_dir / "d1_low_top200_clean_annotated.csv", index=False)

    summary = {
        "annotation_file": str(ann_path),
        "d1_high_summary": summarize_table(d1_high_ann, "d1_high"),
        "d1_low_summary": summarize_table(d1_low_ann, "d1_low"),
    }

    with open(out_dir / "marker_annotation_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))
    print("\n[D1-high annotated top 20]")
    print(d1_high_ann.head(20)[["gene_id", "gene_symbol", "gene_biotype", "delta_mean_d1_minus_d2", "q_value_bh"]])
    print("\n[D1-low annotated top 20]")
    print(d1_low_ann.head(20)[["gene_id", "gene_symbol", "gene_biotype", "delta_mean_d1_minus_d2", "q_value_bh"]])


if __name__ == "__main__":
    main()
