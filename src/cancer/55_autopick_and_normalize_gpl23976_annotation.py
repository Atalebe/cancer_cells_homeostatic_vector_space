#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional

import pandas as pd


def choose_first(cols, hints):
    for hint in hints:
        for c in cols:
            if hint in c.lower():
                return c
    return None


def load_candidate(path: Path):
    text = path.read_text(encoding="utf-8", errors="replace")[:20000]

    if "<table" in text.lower() or "<html" in text.lower():
        tables = pd.read_html(str(path))
        best_df = None
        best_score = -1
        for df in tables:
            cols = [str(c) for c in df.columns]
            score = sum(
                1
                for c in cols
                for h in ["gene", "symbol", "chr", "chrom", "feature", "island", "mapinfo", "id", "probe"]
                if h in c.lower()
            )
            if score > best_score:
                best_score = score
                best_df = df
        if best_df is None:
            raise ValueError("no html table candidate found")
        return best_df.astype(str), "html_table"

    head = text.splitlines()[:20]
    tab_votes = sum("\t" in ln for ln in head)
    comma_votes = sum("," in ln for ln in head)
    sep = "\t" if tab_votes >= comma_votes else ","
    return pd.read_csv(path, sep=sep, dtype=str, low_memory=False), "text_table"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inspection-summary", required=True)
    ap.add_argument("--outcsv", required=True)
    ap.add_argument("--outparquet", required=True)
    ap.add_argument("--summary-json", required=True)
    args = ap.parse_args()

    with open(args.inspection_summary, "r", encoding="utf-8") as f:
        insp = json.load(f)

    candidate_file = insp.get("best_candidate_file")
    if not candidate_file:
        raise SystemExit("no best candidate file found in inspection summary")

    candidate_path = Path(candidate_file)
    df, load_mode = load_candidate(candidate_path)
    cols = [str(c) for c in df.columns]

    probe_col = choose_first(cols, ["probe", "id", "cg", "ilmn", "spot"])
    gene_col = choose_first(cols, ["gene symbol", "gene_symbol", "symbol", "refgene_name", "gene"])
    feature_col = choose_first(cols, ["feature", "relation", "group", "class", "region"])
    island_col = choose_first(cols, ["island"])
    chr_col = choose_first(cols, ["chromosome", "chr"])
    mapinfo_col = choose_first(cols, ["mapinfo", "position", "coord", "start"])
    strand_col = choose_first(cols, ["strand"])
    desc_col = choose_first(cols, ["description", "name", "title"])

    if probe_col is None:
        raise SystemExit("could not determine probe column in best candidate")

    out = pd.DataFrame({"probe_id": df[probe_col].astype(str).str.strip()})

    def add_opt(src: Optional[str], dst: str):
        out[dst] = df[src] if src is not None else pd.NA

    add_opt(gene_col, "gene_symbol")
    add_opt(feature_col, "feature_context")
    add_opt(island_col, "island_context")
    add_opt(chr_col, "chromosome")
    add_opt(mapinfo_col, "mapinfo")
    add_opt(strand_col, "strand")
    add_opt(desc_col, "probe_description")

    out = out.dropna(subset=["probe_id"]).drop_duplicates(subset=["probe_id"]).reset_index(drop=True)

    outcsv = Path(args.outcsv)
    outparquet = Path(args.outparquet)
    summary_json = Path(args.summary_json)
    outcsv.parent.mkdir(parents=True, exist_ok=True)
    outparquet.parent.mkdir(parents=True, exist_ok=True)
    summary_json.parent.mkdir(parents=True, exist_ok=True)

    out.to_csv(outcsv, index=False)
    out.to_parquet(outparquet, index=False)

    summary = {
        "candidate_file": str(candidate_path),
        "load_mode": load_mode,
        "input_rows": int(df.shape[0]),
        "input_columns": cols,
        "normalized_rows": int(out.shape[0]),
        "probe_col_used": probe_col,
        "gene_col_used": gene_col,
        "feature_col_used": feature_col,
        "island_col_used": island_col,
        "chr_col_used": chr_col,
        "mapinfo_col_used": mapinfo_col,
        "strand_col_used": strand_col,
        "desc_col_used": desc_col,
        "n_non_null_gene_symbol": int(out["gene_symbol"].notna().sum()),
        "n_non_null_feature_context": int(out["feature_context"].notna().sum()),
        "n_non_null_island_context": int(out["island_context"].notna().sum()),
        "n_non_null_chromosome": int(out["chromosome"].notna().sum()),
    }
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"[ok] wrote normalized annotation csv: {outcsv}")
    print(f"[ok] wrote normalized annotation parquet: {outparquet}")
    print(f"[ok] wrote summary: {summary_json}")
    print(f"[info] candidate_file: {candidate_path}")
    print(f"[info] probe_col_used: {probe_col}")
    print(f"[info] gene_col_used: {gene_col}")
    print(f"[info] n_non_null_gene_symbol: {summary['n_non_null_gene_symbol']}")


if __name__ == "__main__":
    main()
