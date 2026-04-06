#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

def choose_first(cols: List[str], candidates: List[str]) -> Optional[str]:
    low_map = {c.lower(): c for c in cols}
    for cand in candidates:
        for c in cols:
            if cand in c.lower():
                return c
    return None

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--outcsv", required=True)
    ap.add_argument("--outparquet", required=True)
    ap.add_argument("--summary-json", required=True)
    args = ap.parse_args()

    path = Path(args.input)
    outcsv = Path(args.outcsv)
    outparquet = Path(args.outparquet)
    summary_json = Path(args.summary_json)

    sep = "\t" if path.suffix.lower() in {".tsv", ".txt", ".soft"} else ","
    df = pd.read_csv(path, sep=sep, dtype=str, low_memory=False)

    cols = list(df.columns.astype(str))

    probe_col = choose_first(cols, ["probe", "ilmn", "cg", "id"])
    gene_col = choose_first(cols, ["gene symbol", "gene_symbol", "symbol", "refgene_name", "gene"])
    feature_col = choose_first(cols, ["feature", "relation", "group", "class", "regulatory"])
    island_col = choose_first(cols, ["island"])
    chr_col = choose_first(cols, ["chromosome", "chr"])
    mapinfo_col = choose_first(cols, ["mapinfo", "position", "coord", "start"])
    strand_col = choose_first(cols, ["strand"])
    name_col = choose_first(cols, ["name", "description"])

    if probe_col is None:
        raise SystemExit("could not determine probe column")

    out = pd.DataFrame({
        "probe_id": df[probe_col].astype(str).str.strip()
    })

    def add_optional(src: Optional[str], dst: str) -> None:
        if src is not None:
            out[dst] = df[src]
        else:
            out[dst] = pd.NA

    add_optional(gene_col, "gene_symbol")
    add_optional(feature_col, "feature_context")
    add_optional(island_col, "island_context")
    add_optional(chr_col, "chromosome")
    add_optional(mapinfo_col, "mapinfo")
    add_optional(strand_col, "strand")
    add_optional(name_col, "probe_description")

    out = out.dropna(subset=["probe_id"]).drop_duplicates(subset=["probe_id"]).reset_index(drop=True)

    outcsv.parent.mkdir(parents=True, exist_ok=True)
    outparquet.parent.mkdir(parents=True, exist_ok=True)
    summary_json.parent.mkdir(parents=True, exist_ok=True)

    out.to_csv(outcsv, index=False)
    out.to_parquet(outparquet, index=False)

    summary = {
        "input": str(path),
        "input_rows": int(df.shape[0]),
        "input_cols": cols,
        "normalized_rows": int(out.shape[0]),
        "normalized_cols": list(out.columns),
        "probe_col_used": probe_col,
        "gene_col_used": gene_col,
        "feature_col_used": feature_col,
        "island_col_used": island_col,
        "chr_col_used": chr_col,
        "mapinfo_col_used": mapinfo_col,
        "strand_col_used": strand_col,
        "name_col_used": name_col,
        "n_non_null_gene_symbol": int(out["gene_symbol"].notna().sum()),
        "n_non_null_feature_context": int(out["feature_context"].notna().sum()),
        "n_non_null_island_context": int(out["island_context"].notna().sum()),
        "n_non_null_chromosome": int(out["chromosome"].notna().sum()),
    }
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"[ok] wrote normalized annotation csv: {outcsv}")
    print(f"[ok] wrote normalized annotation parquet: {outparquet}")
    print(f"[ok] wrote summary: {summary_json}")
    print(f"[info] probe_col_used: {probe_col}")
    print(f"[info] gene_col_used: {gene_col}")
    print(f"[info] n_non_null_gene_symbol: {summary['n_non_null_gene_symbol']}")

if __name__ == "__main__":
    main()
