#!/usr/bin/env python3

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
import yaml


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"
ANNOT_PATH = REPO_ROOT / "data" / "reference" / "annotations" / "ensembl_gene_annotation.tsv"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def strip_version(gene_id: str) -> str:
    return str(gene_id).split(".", 1)[0]


def classify_symbol(symbol: str) -> dict[str, bool]:
    s = "" if pd.isna(symbol) else str(symbol)
    return {
        "is_ribo": bool(re.match(r"^(RPL|RPS|MRPL|MRPS)", s)),
        "is_hla": bool(re.match(r"^HLA-", s)),
    }


def annotate_table(df: pd.DataFrame, ann: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["gene_id_stripped"] = out["gene_id"].astype(str).map(strip_version)

    ann2 = ann.copy()
    ann2["gene_id_stripped"] = ann2["gene_id"].astype(str).map(strip_version)

    merged = out.merge(
        ann2[["gene_id_stripped", "gene_symbol", "gene_biotype", "chromosome", "is_mito"]],
        on="gene_id_stripped",
        how="left",
    )

    merged["is_mito"] = merged["is_mito"].fillna(False).astype(bool)

    flags = merged["gene_symbol"].map(classify_symbol).apply(pd.Series)
    merged = pd.concat([merged, flags], axis=1)
    merged["is_ribo"] = merged["is_ribo"].fillna(False).astype(bool)
    merged["is_hla"] = merged["is_hla"].fillna(False).astype(bool)

    return merged


def summarize(df: pd.DataFrame, label: str) -> dict:
    return {
        "label": label,
        "n_rows": int(len(df)),
        "n_with_symbol": int(df["gene_symbol"].notna().sum()),
        "n_mito": int(df["is_mito"].sum()),
        "n_ribo": int(df["is_ribo"].sum()),
        "n_hla": int(df["is_hla"].sum()),
        "top_gene_symbols_first20": df["gene_symbol"].fillna(df["gene_id"]).head(20).tolist(),
    }


def main() -> None:
    cfg = read_config()
    in_dir = REPO_ROOT / cfg["results_dir"] / "d2_marker_scan"
    out_dir = REPO_ROOT / cfg["results_dir"] / "d2_marker_annotation"
    out_dir.mkdir(parents=True, exist_ok=True)

    ann = pd.read_csv(ANNOT_PATH, sep="\t")
    required = {"gene_id", "gene_symbol", "gene_biotype", "chromosome", "is_mito"}
    missing = required - set(ann.columns)
    if missing:
        raise ValueError(f"Annotation file missing columns: {sorted(missing)}")

    high = pd.read_csv(in_dir / "d2_1_high_top200_clean.csv")
    low = pd.read_csv(in_dir / "d2_1_low_top200_clean.csv")

    high_ann = annotate_table(high, ann)
    low_ann = annotate_table(low, ann)

    high_out = out_dir / "d2_1_high_top200_clean_annotated.csv"
    low_out = out_dir / "d2_1_low_top200_clean_annotated.csv"

    high_ann.to_csv(high_out, index=False)
    low_ann.to_csv(low_out, index=False)

    summary = {
        "annotation_file": str(ANNOT_PATH),
        "high_summary": summarize(high_ann, "d2_1_high"),
        "low_summary": summarize(low_ann, "d2_1_low"),
        "outputs": {
            "high_csv": str(high_out),
            "low_csv": str(low_out),
        },
    }

    with open(out_dir / "d2_marker_annotation_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))
    print("\n[D2_1-high annotated top 20]")
    print(high_ann.head(20)[["gene_id", "gene_symbol", "gene_biotype", "delta_mean_a_minus_b", "q_value_bh"]])
    print("\n[D2_1-low annotated top 20]")
    print(low_ann.head(20)[["gene_id", "gene_symbol", "gene_biotype", "delta_mean_a_minus_b", "q_value_bh"]])


if __name__ == "__main__":
    main()
