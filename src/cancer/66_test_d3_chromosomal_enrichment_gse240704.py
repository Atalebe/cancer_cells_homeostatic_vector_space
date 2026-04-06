#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact


CHR_CANDIDATES = [
    "chromosome",
    "CHR",
    "chr",
]

DIRECTION_CANDIDATES = [
    "direction",
]

PROBE_CANDIDATES = [
    "ID_REF",
    "probe_id",
    "ID",
]


def choose_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def clean_chr(x: object) -> str | None:
    if pd.isna(x):
        return None
    s = str(x).strip()
    if not s:
        return None
    s = s.replace("chr", "").replace("CHR", "").strip()
    return s


def fisher_for_membership(
    df: pd.DataFrame,
    chromosome_col: str,
    flag_col: str,
) -> pd.DataFrame:
    rows: list[dict] = []
    total_flag = int(df[flag_col].sum())
    total_not = int((~df[flag_col]).sum())

    for chrom in sorted(df[chromosome_col].dropna().astype(str).unique(), key=lambda x: (len(x), x)):
        in_chr = df[chromosome_col].astype(str) == chrom
        a = int((in_chr & df[flag_col]).sum())
        b = int((in_chr & ~df[flag_col]).sum())
        c = int((~in_chr & df[flag_col]).sum())
        d = int((~in_chr & ~df[flag_col]).sum())

        odds, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        rows.append(
            {
                "chromosome": chrom,
                "flag_n": total_flag,
                "background_n": total_not,
                "in_chr_and_flag": a,
                "in_chr_and_background": b,
                "outside_chr_and_flag": c,
                "outside_chr_and_background": d,
                "frac_flag_in_chr": a / total_flag if total_flag > 0 else None,
                "frac_background_in_chr": b / total_not if total_not > 0 else None,
                "fisher_odds_ratio": odds,
                "fisher_p_value": p,
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["rank_by_p_value"] = out["fisher_p_value"].rank(method="min")
        out = out.sort_values(["fisher_p_value", "chromosome"]).reset_index(drop=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-csv", required=True, help="Directional D3 table, usually D3_vs_not_D3_top_200_probes_with_direction.csv")
    ap.add_argument("--background-csv", required=False, help="Optional background annotation table of all selected probes or all platform probes")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--higher-label", default="higher_in_a")
    ap.add_argument("--lower-label", default="lower_in_a")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input_csv)
    direction_col = choose_col(df, DIRECTION_CANDIDATES)
    chr_col = choose_col(df, CHR_CANDIDATES)
    probe_col = choose_col(df, PROBE_CANDIDATES)

    if direction_col is None:
        raise SystemExit("missing direction column")
    if chr_col is None:
        raise SystemExit("missing chromosome column")
    if probe_col is None:
        raise SystemExit("missing probe column")

    df = df.copy()
    df["chromosome_clean"] = df[chr_col].map(clean_chr)
    df = df[df["chromosome_clean"].notna()].copy()

    higher = df.copy()
    higher["is_target"] = higher[direction_col].astype(str) == args.higher_label
    higher_res = fisher_for_membership(higher, "chromosome_clean", "is_target")
    higher_res.to_csv(outdir / "D3_higher_in_a_chromosome_enrichment.csv", index=False)

    lower = df.copy()
    lower["is_target"] = lower[direction_col].astype(str) == args.lower_label
    lower_res = fisher_for_membership(lower, "chromosome_clean", "is_target")
    lower_res.to_csv(outdir / "D3_lower_in_a_chromosome_enrichment.csv", index=False)

    counts = (
        df.groupby(["chromosome_clean", direction_col], as_index=False)
        .size()
        .rename(columns={"size": "n"})
        .sort_values(["chromosome_clean", direction_col])
        .reset_index(drop=True)
    )
    counts.to_csv(outdir / "D3_directional_chromosome_counts.csv", index=False)

    summary = {
        "input_csv": args.input_csv,
        "direction_col_used": direction_col,
        "chromosome_col_used": chr_col,
        "probe_col_used": probe_col,
        "n_rows_used": int(len(df)),
        "directions_present": sorted(df[direction_col].astype(str).unique().tolist()),
        "n_unique_chromosomes": int(df["chromosome_clean"].nunique()),
        "note": "This is chromosome-level enrichment. Chromosomal-arm enrichment will require cytoband or centromere breakpoint mapping.",
    }
    with open(outdir / "D3_chromosome_enrichment_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"[ok] wrote chromosome enrichment outputs to {outdir}")
    print("[info] summary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
