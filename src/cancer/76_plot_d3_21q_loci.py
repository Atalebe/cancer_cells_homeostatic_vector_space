#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
DATASET = "gse240704"

INDIR = REPO_ROOT / "results" / DATASET / "chr21_locus_followup"
OUTDIR = INDIR / "plots"

SUMMARY_JSON = INDIR / "D3_21q_summary.json"


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def choose_effect_col(df: pd.DataFrame) -> str | None:
    candidates = [
        "delta_beta",
        "mean_diff",
        "effect_size",
        "loading",
        "score",
        "abs_loading",
        "t_stat",
        "logfc",
    ]
    for c in candidates:
        if c in df.columns:
            return c
    return None


def plot_locus_window(csv_path: Path, out_path: Path) -> None:
    df = pd.read_csv(csv_path)
    if df.empty:
        return

    if "position_std" not in df.columns:
        return

    target_gene = str(df["target_gene"].iloc[0]) if "target_gene" in df.columns else "target"
    locus_start = df["window_start"].iloc[0] if "window_start" in df.columns else None
    locus_end = df["window_end"].iloc[0] if "window_end" in df.columns else None

    effect_col = choose_effect_col(df)
    if effect_col is None:
        y = np.ones(len(df))
        ylabel = "probe count"
    else:
        y = pd.to_numeric(df[effect_col], errors="coerce").fillna(0.0)
        ylabel = effect_col

    marker_size = 40
    has_direct = "contains_target_gene_token" in df.columns

    plt.figure(figsize=(12, 4.5))

    if has_direct:
        mask = df["contains_target_gene_token"].fillna(False).astype(bool)
        plt.scatter(
            df.loc[~mask, "position_std"],
            y.loc[~mask] if hasattr(y, "loc") else y[~mask.values],
            alpha=0.8,
            label="other 21q probes",
        )
        plt.scatter(
            df.loc[mask, "position_std"],
            y.loc[mask] if hasattr(y, "loc") else y[mask.values],
            alpha=0.95,
            label=f"{target_gene} annotated probes",
            marker="D",
        )
    else:
        plt.scatter(df["position_std"], y, alpha=0.85)

    if "window_start" in df.columns and "window_end" in df.columns:
        plt.axvspan(
            df["window_start"].iloc[0] + 2_000_000,
            df["window_end"].iloc[0] - 2_000_000,
            alpha=0.12,
        )

    if "distance_to_locus_start_bp" in df.columns:
        pass

    title = f"D3 lower-in-a 21q probe window around {target_gene}"
    plt.title(title)
    plt.xlabel("Genomic position on chr21 (hg19)")
    plt.ylabel(ylabel)

    if has_direct and df["contains_target_gene_token"].fillna(False).any():
        plt.legend(frameon=False)

    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def main() -> None:
    ensure_dir(OUTDIR)

    if not SUMMARY_JSON.exists():
        raise FileNotFoundError(f"Missing summary file: {SUMMARY_JSON}")

    with open(SUMMARY_JSON, "r", encoding="utf-8") as f:
        summary = json.load(f)

    locus_windows = summary.get("outputs", {}).get("locus_windows", {})
    for gene, csv_path in locus_windows.items():
        csv_path = Path(csv_path)
        if not csv_path.exists():
            continue
        out_path = OUTDIR / f"D3_21q_{gene}_probe_window.png"
        plot_locus_window(csv_path=csv_path, out_path=out_path)

    print(f"[ok] wrote plots to {OUTDIR}")


if __name__ == "__main__":
    main()
