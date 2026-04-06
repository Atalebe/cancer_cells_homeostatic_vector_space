#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def safe_log2_or(x):
    if pd.isna(x):
        return np.nan
    if x == 0:
        return -6.0
    if np.isinf(x):
        return 6.0
    if x < 0:
        return np.nan
    return np.log2(x)


def chr_sort_key(ch):
    s = str(ch)
    if s == "X":
        return (1000, s)
    if s == "Y":
        return (1001, s)
    try:
        return (int(s), s)
    except Exception:
        return (999, s)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--enrichment-csv", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.enrichment_csv)
    required = {"direction", "chromosome", "fisher_odds_ratio", "fisher_p_value"}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"missing required columns: {sorted(missing)}")

    df = df.copy()
    df["log2_or_capped"] = df["fisher_odds_ratio"].map(safe_log2_or)
    df["neglog10_p"] = -np.log10(df["fisher_p_value"].clip(lower=1e-300))

    chroms = sorted(df["chromosome"].dropna().astype(str).unique().tolist(), key=chr_sort_key)
    directions = sorted(df["direction"].dropna().astype(str).unique().tolist())

    pivot_or = (
        df.assign(chromosome=df["chromosome"].astype(str))
        .pivot(index="chromosome", columns="direction", values="log2_or_capped")
        .reindex(chroms)
    )
    pivot_p = (
        df.assign(chromosome=df["chromosome"].astype(str))
        .pivot(index="chromosome", columns="direction", values="neglog10_p")
        .reindex(chroms)
    )

    # Plot 1: directional log2 odds ratios
    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(chroms))
    width = 0.38 if len(directions) >= 2 else 0.6

    offsets = np.linspace(-width / 2, width / 2, max(len(directions), 1))
    for i, direction in enumerate(directions):
        vals = pivot_or[direction].values if direction in pivot_or.columns else np.full(len(chroms), np.nan)
        ax.bar(x + offsets[i], vals, width=width / max(len(directions), 1), label=direction)

    ax.axhline(0, linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(chroms, rotation=0)
    ax.set_ylabel("log2 Fisher odds ratio, capped")
    ax.set_xlabel("Chromosome")
    ax.set_title("D3 directional chromosome enrichment vs selected-probe background")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "D3_directional_chromosome_log2OR_vs_selected_background.png", dpi=200)
    plt.close(fig)

    # Plot 2: significance
    fig, ax = plt.subplots(figsize=(12, 5))
    for i, direction in enumerate(directions):
        vals = pivot_p[direction].values if direction in pivot_p.columns else np.full(len(chroms), np.nan)
        ax.bar(x + offsets[i], vals, width=width / max(len(directions), 1), label=direction)

    ax.axhline(-np.log10(0.05), linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(chroms, rotation=0)
    ax.set_ylabel("-log10 p")
    ax.set_xlabel("Chromosome")
    ax.set_title("D3 directional chromosome enrichment significance")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "D3_directional_chromosome_neglog10p_vs_selected_background.png", dpi=200)
    plt.close(fig)

    # Plot 3: chromosome 21 zoom
    df21 = df.loc[df["chromosome"].astype(str) == "21"].copy()
    if not df21.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        x2 = np.arange(df21.shape[0])
        ax.bar(x2, df21["log2_or_capped"].values)
        ax.set_xticks(x2)
        ax.set_xticklabels(df21["direction"].tolist(), rotation=15)
        ax.axhline(0, linewidth=1)
        ax.set_ylabel("log2 Fisher odds ratio, capped")
        ax.set_title("Chromosome 21 directional enrichment")
        fig.tight_layout()
        fig.savefig(outdir / "D3_chromosome21_directional_log2OR.png", dpi=200)
        plt.close(fig)

    print("[ok] wrote enrichment plots to", outdir)
    print("[info] files:")
    print("  D3_directional_chromosome_log2OR_vs_selected_background.png")
    print("  D3_directional_chromosome_neglog10p_vs_selected_background.png")
    if not df21.empty:
        print("  D3_chromosome21_directional_log2OR.png")


if __name__ == "__main__":
    main()
