#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot D3 diagnostics.")
    p.add_argument("--sample-metrics-csv", required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def boxplot(df: pd.DataFrame, xcol: str, ycol: str, outpath: Path, title: str) -> None:
    sub = df[[xcol, ycol]].copy()
    sub[ycol] = pd.to_numeric(sub[ycol], errors="coerce")
    sub = sub.dropna()
    groups = []
    labels = []
    for label in ["D1", "D2", "D3"]:
        vals = sub.loc[sub[xcol] == label, ycol].tolist()
        if vals:
            groups.append(vals)
            labels.append(label)
    if not groups:
        return
    plt.figure(figsize=(6, 4.5))
    plt.boxplot(groups, labels=labels)
    plt.title(title)
    plt.xlabel(xcol)
    plt.ylabel(ycol)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def scatter(df: pd.DataFrame, x: str, y: str, outpath: Path, title: str) -> None:
    sub = df[[x, y, "state_domain"]].copy()
    sub[x] = pd.to_numeric(sub[x], errors="coerce")
    sub[y] = pd.to_numeric(sub[y], errors="coerce")
    sub = sub.dropna()
    plt.figure(figsize=(6, 5))
    for dom in ["D1", "D2", "D3"]:
        d = sub[sub["state_domain"] == dom]
        if len(d):
            plt.scatter(d[x], d[y], label=dom, alpha=0.7, s=20)
    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.sample_metrics_csv)

    boxplot(df, "state_domain", "mean_beta", outdir / "mean_beta_by_state_domain.png", "Mean beta by state domain")
    boxplot(df, "state_domain", "std_beta", outdir / "std_beta_by_state_domain.png", "Beta standard deviation by state domain")
    boxplot(df, "state_domain", "iqr_beta", outdir / "iqr_beta_by_state_domain.png", "Beta IQR by state domain")
    boxplot(df, "state_domain", "mean_neighbor_distance", outdir / "mean_neighbor_distance_by_state_domain.png", "Mean neighbor distance by state domain")
    boxplot(df, "state_domain", "min_hrsm_distance_to_shell", outdir / "min_hrsm_distance_to_shell_by_state_domain.png", "Min HRSM distance to shell by state domain")
    scatter(df, "mean_beta", "std_beta", outdir / "mean_beta_vs_std_beta_by_domain.png", "Mean beta vs beta standard deviation")
    scatter(df, "phi", "std_beta", outdir / "phi_vs_std_beta_by_domain.png", "phi vs beta standard deviation")
    scatter(df, "S", "std_beta", outdir / "S_vs_std_beta_by_domain.png", "S vs beta standard deviation")
    scatter(df, "phi", "min_hrsm_distance_to_shell", outdir / "phi_vs_min_hrsm_distance_to_shell.png", "phi vs min shell distance")

    print("[ok] wrote D3 diagnostic plots to", outdir)


if __name__ == "__main__":
    main()
