#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def choose_probe_col(df: pd.DataFrame) -> str | None:
    preferred = [
        "ID",
        "ID_REF",
        "Name",
        "IlmnID",
        "Probe_ID",
        "probe_id",
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state-table", required=True)
    ap.add_argument("--pca-scores-csv", required=True)
    ap.add_argument("--top-probe-loadings-csv", required=True)
    ap.add_argument("--annotation-csv", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--top-n", type=int, default=200)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    state = pd.read_parquet(args.state_table)
    scores = pd.read_csv(args.pca_scores_csv)
    loadings = pd.read_csv(args.top_probe_loadings_csv)
    ann = pd.read_csv(args.annotation_csv, dtype=str, low_memory=False)

    merged = state.merge(scores, on="sample_id", how="inner")
    pc_cols = [c for c in scores.columns if c.startswith("PC")]

    rows = []
    for pc in pc_cols:
        x = pd.to_numeric(merged[pc], errors="coerce")
        y = pd.to_numeric(merged["S"], errors="coerce")
        mask = x.notna() & y.notna()
        if mask.sum() < 5:
            continue
        corr = np.corrcoef(x[mask], y[mask])[0, 1]
        rows.append(
            {
                "pc": pc,
                "corr_with_S": float(corr),
                "abs_corr_with_S": float(abs(corr)),
                "n": int(mask.sum()),
            }
        )

    corr_df = pd.DataFrame(rows).sort_values("abs_corr_with_S", ascending=False).reset_index(drop=True)
    corr_df.to_csv(outdir / "s_axis_pc_correlation_table.csv", index=False)

    if corr_df.empty:
        raise SystemExit("no usable PC correlations with S were found")

    best_pc = corr_df.iloc[0]["pc"]
    best_pc_label = best_pc.upper()

    sub = loadings[loadings["pc"].str.upper() == best_pc_label].copy()
    if sub.empty:
        raise SystemExit(f"no loading rows found for {best_pc_label}")

    sub = sub.sort_values("abs_loading", ascending=False).head(args.top_n).copy()

    probe_col_ann = choose_probe_col(ann)
    if probe_col_ann is None:
        raise SystemExit("could not determine probe column in parsed annotation table")

    keep_cols = [probe_col_ann]
    for c in [
        "Gene_Symbol",
        "UCSC_RefGene_Name",
        "UCSC_RefGene_Group",
        "Relation_to_Island",
        "CHR",
        "MAPINFO",
        "Chromosome",
        "Coordinate",
        "Description",
        "Gene",
        "Symbol",
    ]:
        if c in ann.columns and c not in keep_cols:
            keep_cols.append(c)

    ann2 = ann[keep_cols].drop_duplicates(subset=[probe_col_ann]).copy()
    ann2[probe_col_ann] = ann2[probe_col_ann].astype(str)
    sub["ID_REF"] = sub["ID_REF"].astype(str)

    merged_top = sub.merge(ann2, left_on="ID_REF", right_on=probe_col_ann, how="left")
    merged_top.to_csv(outdir / "s_axis_top_probe_drivers_annotated.csv", index=False)

    summary = {
        "best_pc_for_S": best_pc,
        "corr_with_S": float(corr_df.iloc[0]["corr_with_S"]),
        "abs_corr_with_S": float(corr_df.iloc[0]["abs_corr_with_S"]),
        "top_n": int(args.top_n),
        "annotation_probe_col": probe_col_ann,
        "n_annotated_rows": int(merged_top[probe_col_ann].notna().sum()),
    }
    with open(outdir / "s_axis_driver_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote S-axis PC correlation table:", outdir / "s_axis_pc_correlation_table.csv")
    print("[ok] wrote annotated S-axis top probes:", outdir / "s_axis_top_probe_drivers_annotated.csv")
    print("[ok] wrote summary:", outdir / "s_axis_driver_summary.json")
    print("[info] best_pc_for_S:", best_pc)
    print("[info] corr_with_S:", corr_df.iloc[0]["corr_with_S"])


if __name__ == "__main__":
    main()
