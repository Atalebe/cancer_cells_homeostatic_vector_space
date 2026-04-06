#!/usr/bin/env python3

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
GTF_PATH = REPO_ROOT / "data" / "reference" / "gencode" / "gencode.v49.basic.annotation.gtf.gz"
OUTDIR = REPO_ROOT / "data" / "reference" / "annotations"
OUTPATH = OUTDIR / "ensembl_gene_annotation.tsv"


def parse_attrs(attr_text: str) -> dict[str, str]:
    out = {}
    for key, val in re.findall(r'(\S+)\s+"([^"]+)"', attr_text):
        out[key] = val
    return out


def main() -> None:
    if not GTF_PATH.exists():
        raise FileNotFoundError(f"Missing GTF: {GTF_PATH}")

    OUTDIR.mkdir(parents=True, exist_ok=True)

    rows = []
    with gzip.open(GTF_PATH, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature != "gene":
                continue

            d = parse_attrs(attrs)
            gene_id = d.get("gene_id")
            gene_symbol = d.get("gene_name")
            gene_biotype = d.get("gene_type") or d.get("gene_biotype")

            if not gene_id:
                continue

            rows.append(
                {
                    "gene_id": gene_id,
                    "gene_symbol": gene_symbol,
                    "gene_biotype": gene_biotype,
                    "chromosome": chrom,
                    "is_mito": chrom in {"chrM", "MT"},
                }
            )

    df = pd.DataFrame(rows).drop_duplicates(subset=["gene_id"]).reset_index(drop=True)
    df.to_csv(OUTPATH, sep="\t", index=False)

    print(f"[ok] wrote {OUTPATH}")
    print(df.head())


if __name__ == "__main__":
    main()
