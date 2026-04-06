from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def pick_col(columns, candidates):
    lower_map = {c.lower(): c for c in columns}
    for cand in candidates:
        for c in columns:
            if cand in c.lower():
                return c
    return None


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare a curated metadata mapping template from series matrix sample metadata."
    )
    parser.add_argument("--series-metadata-csv", required=True)
    parser.add_argument("--mapping-template-csv", required=True)
    parser.add_argument("--outcsv", required=True)
    args = parser.parse_args()

    sm = pd.read_csv(args.series_metadata_csv, dtype=str).fillna("")
    mt = pd.read_csv(args.mapping_template_csv, dtype=str).fillna("")

    if "sample_id" not in sm.columns or "sample_id" not in mt.columns:
        raise ValueError("Both inputs must contain sample_id")

    merged = mt.merge(sm, on="sample_id", how="left", suffixes=("", "_series"))

    title_col = pick_col(sm.columns, ["title"])
    source_col = pick_col(sm.columns, ["source_name"])
    char_cols = [c for c in sm.columns if "characteristics" in c.lower()]
    desc_col = pick_col(sm.columns, ["description"])
    supp_col = pick_col(sm.columns, ["supplementary_file"])

    if title_col and "biological_condition" not in merged.columns:
        merged["biological_condition"] = merged[title_col]
    if source_col and "cell_type" not in merged.columns:
        merged["cell_type"] = merged[source_col]
    if desc_col and "source_note" not in merged.columns:
        merged["source_note"] = merged[desc_col]
    if supp_col and "supplementary_file" not in merged.columns:
        merged["supplementary_file"] = merged[supp_col]

    for i, c in enumerate(char_cols, start=1):
        merged[f"characteristics_ch{i}"] = merged[c]

    Path(args.outcsv).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.outcsv, index=False)

    print("[ok] wrote prepared metadata mapping template:", args.outcsv)
    print("[info] rows:", len(merged))
    print("[info] title_col:", title_col)
    print("[info] source_col:", source_col)
    print("[info] desc_col:", desc_col)
    print("[info] characteristic columns:", len(char_cols))


if __name__ == "__main__":
    main()
