#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def find_col(df: pd.DataFrame, candidates: list[str]) -> Optional[str]:
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    return None


def clean_chr(val) -> Optional[str]:
    if pd.isna(val):
        return None
    s = str(val).strip()
    if not s:
        return None
    s = s.replace("chr", "").replace("CHR", "").strip()
    if s == "23":
        return "X"
    if s == "24":
        return "Y"
    if s.upper() in {"M", "MT"}:
        return None
    return s


def assign_arm(chromosome: str, mapinfo, centromere_df: pd.DataFrame) -> Optional[str]:
    if chromosome is None or pd.isna(mapinfo):
        return None
    try:
        pos = float(mapinfo)
    except Exception:
        return None

    row = centromere_df.loc[centromere_df["chromosome"].astype(str) == str(chromosome)]
    if row.empty:
        return None
    row = row.iloc[0]

    p_end = float(row["p_end"])
    q_start = float(row["q_start"])

    if pos <= p_end:
        return f"{chromosome}p"
    if pos >= q_start:
        return f"{chromosome}q"
    return f"{chromosome}cen"


def compute_enrichment(flag_arm: pd.Series, bg_arm: pd.Series, direction_label: str) -> pd.DataFrame:
    flag_arm = flag_arm.dropna().astype(str)
    bg_arm = bg_arm.dropna().astype(str)

    arms = sorted(set(flag_arm.unique()) | set(bg_arm.unique()), key=str)

    rows = []
    flag_n = int(flag_arm.shape[0])
    bg_n = int(bg_arm.shape[0])

    for arm in arms:
        a = int((flag_arm == arm).sum())
        b = int((bg_arm == arm).sum())
        c = flag_n - a
        d = bg_n - b

        odds_ratio, p_value = fisher_exact([[a, c], [b, d]], alternative="two-sided")

        rows.append(
            {
                "direction": direction_label,
                "chromosome_arm": arm,
                "flag_n": flag_n,
                "background_n": bg_n,
                "in_arm_and_flag": a,
                "in_arm_and_background": b,
                "outside_arm_and_flag": c,
                "outside_arm_and_background": d,
                "frac_flag_in_arm": a / flag_n if flag_n else np.nan,
                "frac_background_in_arm": b / bg_n if bg_n else np.nan,
                "fisher_odds_ratio": odds_ratio,
                "fisher_p_value": p_value,
            }
        )

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["fisher_p_value", "chromosome_arm"]).reset_index(drop=True)
        out["rank_by_p_value"] = np.arange(1, len(out) + 1)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--directional-csv", required=True)
    ap.add_argument("--background-annotation-parquet", required=True)
    ap.add_argument("--selected-probe-matrix", required=True)
    ap.add_argument("--centromere-csv", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    d3 = pd.read_csv(args.directional_csv)
    ann = pd.read_parquet(args.background_annotation_parquet)
    selected = pd.read_parquet(args.selected_probe_matrix)
    cent = pd.read_csv(args.centromere_csv)

    needed_cent_cols = {"chromosome", "p_end", "q_start"}
    missing = needed_cent_cols - set(cent.columns)
    if missing:
        raise RuntimeError(f"centromere csv missing required columns: {sorted(missing)}")

    probe_col_ann = find_col(ann, ["ID_REF", "probe_id", "ID"])
    chr_col_ann = find_col(ann, ["chromosome", "CHR", "chr"])
    map_col_ann = find_col(ann, ["mapinfo", "MAPINFO"])

    if probe_col_ann is None or chr_col_ann is None or map_col_ann is None:
        raise RuntimeError("could not determine probe/chromosome/mapinfo columns in annotation parquet")

    probe_col_selected = find_col(selected, ["ID_REF", "probe_id", "ID"])
    if probe_col_selected is None:
        raise RuntimeError("could not determine probe column in selected probe matrix")

    probe_col_d3 = find_col(d3, ["ID_REF", "probe_id", "ID"])
    dir_col_d3 = find_col(d3, ["direction"])
    if probe_col_d3 is None or dir_col_d3 is None:
        raise RuntimeError("could not determine probe or direction columns in directional csv")

    ann_sub = ann[[probe_col_ann, chr_col_ann, map_col_ann]].copy()
    ann_sub = ann_sub.rename(
        columns={
            probe_col_ann: "ID_REF",
            chr_col_ann: "chromosome",
            map_col_ann: "mapinfo",
        }
    )
    ann_sub["chromosome"] = ann_sub["chromosome"].map(clean_chr)
    ann_sub["chromosome_arm"] = ann_sub.apply(
        lambda r: assign_arm(r["chromosome"], r["mapinfo"], cent), axis=1
    )
    ann_sub = ann_sub.dropna(subset=["ID_REF"]).drop_duplicates(subset=["ID_REF"])

    selected_sub = selected[[probe_col_selected]].copy().rename(columns={probe_col_selected: "ID_REF"})
    selected_sub = selected_sub.dropna(subset=["ID_REF"]).drop_duplicates(subset=["ID_REF"])

    bg = selected_sub.merge(
        ann_sub[["ID_REF", "chromosome", "mapinfo", "chromosome_arm"]],
        on="ID_REF",
        how="left",
    )
    bg = bg.dropna(subset=["chromosome_arm"]).copy()

    d3_sub = d3[[probe_col_d3, dir_col_d3]].copy().rename(
        columns={probe_col_d3: "ID_REF", dir_col_d3: "direction"}
    )
    d3_sub = d3_sub.dropna(subset=["ID_REF", "direction"]).drop_duplicates(subset=["ID_REF", "direction"])
    d3_sub = d3_sub.merge(
        ann_sub[["ID_REF", "chromosome", "mapinfo", "chromosome_arm"]],
        on="ID_REF",
        how="left",
    )
    d3_sub = d3_sub.dropna(subset=["chromosome_arm"]).copy()

    merged_out = outdir / "D3_directional_probe_arm_assignments.csv"
    d3_sub.to_csv(merged_out, index=False)

    all_rows = []
    for direction in sorted(d3_sub["direction"].dropna().unique()):
        flag = d3_sub.loc[d3_sub["direction"] == direction].copy()
        flag_ids = set(flag["ID_REF"].astype(str))
        bg_rest = bg.loc[~bg["ID_REF"].astype(str).isin(flag_ids)].copy()

        res = compute_enrichment(
            flag_arm=flag["chromosome_arm"],
            bg_arm=bg_rest["chromosome_arm"],
            direction_label=str(direction),
        )
        outfile = outdir / f"D3_{direction}_chromosome_arm_enrichment.csv"
        res.to_csv(outfile, index=False)
        all_rows.append(res)

    combined = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    combined_out = outdir / "D3_directional_chromosome_arm_enrichment_combined.csv"
    combined.to_csv(combined_out, index=False)

    counts = (
        d3_sub.groupby(["direction", "chromosome_arm"])
        .size()
        .reset_index(name="n")
        .sort_values(["direction", "chromosome_arm"])
        .reset_index(drop=True)
    )
    counts_out = outdir / "D3_directional_chromosome_arm_counts.csv"
    counts.to_csv(counts_out, index=False)

    summary = {
        "directional_csv": args.directional_csv,
        "background_annotation_parquet": args.background_annotation_parquet,
        "selected_probe_matrix": args.selected_probe_matrix,
        "centromere_csv": args.centromere_csv,
        "n_d3_rows_with_arm": int(d3_sub.shape[0]),
        "n_selected_background_with_arm": int(bg.shape[0]),
        "directions_present": sorted(d3_sub["direction"].dropna().unique().tolist()),
        "note": "Arm assignments are based on chromosome and MAPINFO against supplied centromere breakpoints.",
    }
    with open(outdir / "D3_directional_chromosome_arm_enrichment_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote chromosome-arm enrichment outputs to", outdir)
    print("[info] summary:")
    print("  n_d3_rows_with_arm:", summary["n_d3_rows_with_arm"])
    print("  n_selected_background_with_arm:", summary["n_selected_background_with_arm"])
    print("  directions_present:", summary["directions_present"])


if __name__ == "__main__":
    main()
