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
    if s in {"23"}:
        return "X"
    if s in {"24"}:
        return "Y"
    if s in {"M", "MT", "m", "mt"}:
        return None
    return s


def compute_one_direction(
    flag_chr_series: pd.Series,
    bg_chr_series: pd.Series,
    direction_label: str,
) -> pd.DataFrame:
    rows = []

    bg_chr_series = bg_chr_series.dropna().astype(str)
    flag_chr_series = flag_chr_series.dropna().astype(str)

    chromosomes = sorted(set(bg_chr_series.unique()) | set(flag_chr_series.unique()), key=str)

    flag_n = int(flag_chr_series.shape[0])
    bg_n = int(bg_chr_series.shape[0])

    for chrom in chromosomes:
        in_chr_and_flag = int((flag_chr_series == chrom).sum())
        in_chr_and_background = int((bg_chr_series == chrom).sum())

        outside_chr_and_flag = flag_n - in_chr_and_flag
        outside_chr_and_background = bg_n - in_chr_and_background

        table = np.array(
            [
                [in_chr_and_flag, outside_chr_and_flag],
                [in_chr_and_background, outside_chr_and_background],
            ],
            dtype=int,
        )

        odds_ratio, p_value = fisher_exact(table, alternative="two-sided")

        frac_flag_in_chr = (
            in_chr_and_flag / flag_n if flag_n > 0 else np.nan
        )
        frac_background_in_chr = (
            in_chr_and_background / bg_n if bg_n > 0 else np.nan
        )

        rows.append(
            {
                "direction": direction_label,
                "chromosome": chrom,
                "flag_n": flag_n,
                "background_n": bg_n,
                "in_chr_and_flag": in_chr_and_flag,
                "in_chr_and_background": in_chr_and_background,
                "outside_chr_and_flag": outside_chr_and_flag,
                "outside_chr_and_background": outside_chr_and_background,
                "frac_flag_in_chr": frac_flag_in_chr,
                "frac_background_in_chr": frac_background_in_chr,
                "fisher_odds_ratio": odds_ratio,
                "fisher_p_value": p_value,
            }
        )

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["fisher_p_value", "chromosome"]).reset_index(drop=True)
        out["rank_by_p_value"] = np.arange(1, len(out) + 1)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--directional-csv", required=True)
    ap.add_argument("--background-annotation-parquet", required=True)
    ap.add_argument("--selected-probe-matrix", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    d3 = pd.read_csv(args.directional_csv)
    ann = pd.read_parquet(args.background_annotation_parquet)

    # Use selected probe universe as the background
    probe_matrix = pd.read_parquet(args.selected_probe_matrix)
    probe_col_matrix = find_col(probe_matrix, ["ID_REF", "probe_id", "ID"])
    if probe_col_matrix is None:
        raise RuntimeError("could not determine probe column in selected probe matrix")

    probe_col_ann = find_col(ann, ["ID_REF", "probe_id", "ID"])
    if probe_col_ann is None:
        raise RuntimeError("could not determine probe column in background annotation parquet")

    chr_col_ann = find_col(ann, ["chromosome", "CHR", "chr"])
    if chr_col_ann is None:
        raise RuntimeError("could not determine chromosome column in background annotation parquet")

    probe_col_d3 = find_col(d3, ["ID_REF", "probe_id", "ID"])
    if probe_col_d3 is None:
        raise RuntimeError("could not determine probe column in directional csv")

    direction_col = find_col(d3, ["direction"])
    if direction_col is None:
        raise RuntimeError("could not determine direction column in directional csv")

    ann_sub = ann[[probe_col_ann, chr_col_ann]].copy()
    ann_sub = ann_sub.rename(
        columns={
            probe_col_ann: "ID_REF",
            chr_col_ann: "chromosome",
        }
    )
    ann_sub["chromosome_clean"] = ann_sub["chromosome"].map(clean_chr)
    ann_sub = ann_sub.dropna(subset=["ID_REF"]).drop_duplicates(subset=["ID_REF"])

    selected = probe_matrix[[probe_col_matrix]].copy()
    selected = selected.rename(columns={probe_col_matrix: "ID_REF"})
    selected = selected.dropna(subset=["ID_REF"]).drop_duplicates(subset=["ID_REF"])

    bg = selected.merge(ann_sub[["ID_REF", "chromosome_clean"]], on="ID_REF", how="left")
    bg = bg.dropna(subset=["chromosome_clean"]).copy()

    d3_sub = d3[[probe_col_d3, direction_col]].copy()
    d3_sub = d3_sub.rename(columns={probe_col_d3: "ID_REF", direction_col: "direction"})
    d3_sub = d3_sub.dropna(subset=["ID_REF", "direction"]).drop_duplicates(subset=["ID_REF", "direction"])
    d3_sub = d3_sub.merge(ann_sub[["ID_REF", "chromosome_clean"]], on="ID_REF", how="left")
    d3_sub = d3_sub.dropna(subset=["chromosome_clean"]).copy()

    # Background for each direction is selected probes excluding that direction's flagged probes
    combined_rows = []
    written_files = []

    for direction in sorted(d3_sub["direction"].dropna().unique()):
        flag = d3_sub.loc[d3_sub["direction"] == direction, ["ID_REF", "chromosome_clean"]].copy()
        flag_probe_ids = set(flag["ID_REF"].astype(str))

        bg_rest = bg.loc[~bg["ID_REF"].astype(str).isin(flag_probe_ids), ["ID_REF", "chromosome_clean"]].copy()

        res = compute_one_direction(
            flag_chr_series=flag["chromosome_clean"],
            bg_chr_series=bg_rest["chromosome_clean"],
            direction_label=str(direction),
        )
        combined_rows.append(res)

        outfile = outdir / f"D3_{direction}_chromosome_enrichment_vs_selected_background.csv"
        res.to_csv(outfile, index=False)
        written_files.append(str(outfile))

    combined = pd.concat(combined_rows, ignore_index=True) if combined_rows else pd.DataFrame()
    combined_out = outdir / "D3_directional_chromosome_enrichment_vs_selected_background.csv"
    combined.to_csv(combined_out, index=False)
    written_files.append(str(combined_out))

    counts = (
        d3_sub.groupby(["chromosome_clean", "direction"])
        .size()
        .reset_index(name="n")
        .sort_values(["chromosome_clean", "direction"])
        .reset_index(drop=True)
    )
    counts_out = outdir / "D3_directional_chromosome_counts_vs_selected_background.csv"
    counts.to_csv(counts_out, index=False)
    written_files.append(str(counts_out))

    summary = {
        "directional_csv": args.directional_csv,
        "background_annotation_parquet": args.background_annotation_parquet,
        "selected_probe_matrix": args.selected_probe_matrix,
        "n_d3_rows_input": int(d3.shape[0]),
        "n_d3_rows_with_chr": int(d3_sub.shape[0]),
        "n_selected_background_probes": int(bg.shape[0]),
        "direction_col_used": "direction",
        "probe_col_used_d3": "ID_REF",
        "probe_col_used_background": "ID_REF",
        "chromosome_col_used_background": "chromosome",
        "directions_present": sorted(d3_sub["direction"].dropna().unique().tolist()),
        "files_written": written_files,
        "note": "This rerun uses the selected probe universe as the chromosome background, excluding flagged probes for each direction.",
    }
    with open(outdir / "D3_directional_chromosome_enrichment_vs_selected_background_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote chromosome enrichment vs selected background to", outdir)
    print("[info] summary:")
    print("  n_d3_rows_with_chr:", summary["n_d3_rows_with_chr"])
    print("  n_selected_background_probes:", summary["n_selected_background_probes"])
    print("  directions_present:", summary["directions_present"])


if __name__ == "__main__":
    main()
