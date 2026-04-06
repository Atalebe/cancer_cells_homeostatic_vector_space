#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def find_best_pc_for_s(state: pd.DataFrame):
    pc_cols = [c for c in state.columns if c.upper().startswith("PC")]
    out = []
    if "S" not in state.columns:
        return pd.DataFrame()

    sub = state[["S"] + pc_cols].copy()
    for col in pc_cols:
        x = pd.to_numeric(sub[col], errors="coerce")
        y = pd.to_numeric(sub["S"], errors="coerce")
        mask = x.notna() & y.notna()
        if mask.sum() < 5:
            continue
        corr = np.corrcoef(x[mask], y[mask])[0, 1]
        out.append(
            {
                "pc": col,
                "corr_with_S": float(corr),
                "abs_corr_with_S": float(abs(corr)),
                "n": int(mask.sum()),
            }
        )
    df = pd.DataFrame(out)
    if not df.empty:
        df = df.sort_values("abs_corr_with_S", ascending=False).reset_index(drop=True)
    return df


def choose_probe_col(df: pd.DataFrame) -> str | None:
    for col in ["ID_REF", "probe_id", "Probe_ID", "Name", "IlmnID"]:
        if col in df.columns:
            return col
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state-table", required=True)
    ap.add_argument("--top-probe-loadings", required=True)
    ap.add_argument("--annotation-parquet", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--top-n", type=int, default=200)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    state = pd.read_parquet(args.state_table)
    load = pd.read_csv(args.top_probe_loadings)
    ann = pd.read_parquet(args.annotation_parquet)

    pc_map = find_best_pc_for_s(state)
    pc_map.to_csv(outdir / "s_axis_pc_correlation_table.csv", index=False)

    if pc_map.empty:
        raise SystemExit("no PC columns found for S-axis mapping")

    best_pc = pc_map.iloc[0]["pc"]

    sub = load[load["pc"] == best_pc].copy()
    sub = sub.sort_values("abs_loading", ascending=False).head(args.top_n).copy()

    probe_col_ann = choose_probe_col(ann)
    if probe_col_ann is None:
        raise SystemExit("could not determine probe column in annotation parquet")

    keep_cols = [probe_col_ann]
    for c in [
        "Gene_Symbol",
        "UCSC_RefGene_Name",
        "UCSC_RefGene_Group",
        "Relation_to_Island",
        "CHR",
        "MAPINFO",
        "Description",
        "Name",
    ]:
        if c in ann.columns and c not in keep_cols:
            keep_cols.append(c)

    ann2 = ann[keep_cols].drop_duplicates(subset=[probe_col_ann]).copy()
    ann2[probe_col_ann] = ann2[probe_col_ann].astype(str)
    sub["ID_REF"] = sub["ID_REF"].astype(str)

    merged = sub.merge(ann2, left_on="ID_REF", right_on=probe_col_ann, how="left")
    merged.to_csv(outdir / "s_axis_top_probe_drivers_annotated.csv", index=False)

    summary = {
        "best_pc_for_S": best_pc,
        "corr_with_S": float(pc_map.iloc[0]["corr_with_S"]),
        "abs_corr_with_S": float(pc_map.iloc[0]["abs_corr_with_S"]),
        "top_n": int(args.top_n),
        "annotation_probe_col": probe_col_ann,
    }
    with open(outdir / "s_axis_driver_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote S-axis PC mapping:", outdir / "s_axis_pc_correlation_table.csv")
    print("[ok] wrote S-axis top probes:", outdir / "s_axis_top_probe_drivers_annotated.csv")
    print("[ok] wrote summary:", outdir / "s_axis_driver_summary.json")
    print("[info] best_pc_for_S:", best_pc)
    print("[info] corr_with_S:", pc_map.iloc[0]["corr_with_S"])


if __name__ == "__main__":
    main()
