#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


AXES = ["H", "S", "M", "R", "phi"]


def ensure_outdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def load_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.suffix == ".parquet":
        return pd.read_parquet(p)
    return pd.read_csv(p)


def normalize_condition_value(x: object) -> object:
    if pd.isna(x):
        return None
    s = str(x).strip().lower()
    if s in {"mgmt_methylated", "methylated"}:
        return "mgmt_methylated"
    if s in {"mgmt_non_methylated", "non_methylated", "non-methylated"}:
        return "mgmt_non_methylated"
    if "non" in s and "methyl" in s:
        return "mgmt_non_methylated"
    if "methyl" in s:
        return "mgmt_methylated"
    return None


def choose_condition_column(df: pd.DataFrame) -> str:
    for col in [
        "placeholder_condition_cur",
        "biological_condition",
        "tumor_status",
        "placeholder_condition",
        "placeholder_condition_ann",
    ]:
        if col in df.columns and df[col].notna().sum() > 0:
            return col
    raise ValueError(f"No usable condition column found. Available columns: {list(df.columns)}")


def resolve_axis_column(df: pd.DataFrame, axis: str) -> str:
    for col in [axis, f"{axis}_dom", f"{axis}_ann", f"{axis}_state", f"{axis}_state_tbl"]:
        if col in df.columns:
            return col
    raise KeyError(f"Could not resolve axis column for {axis}. Available columns: {list(df.columns)}")


def boxplot_by_domain(df: pd.DataFrame, axis_col: str, axis_label: str, outpath: Path) -> None:
    order = sorted(df["state_domain"].dropna().unique())
    data = [pd.to_numeric(df.loc[df["state_domain"] == d, axis_col], errors="coerce").dropna() for d in order]

    plt.figure(figsize=(7, 5))
    plt.boxplot(data, tick_labels=order)
    plt.xlabel("State domain")
    plt.ylabel(axis_label)
    plt.title(f"{axis_label} by state domain")
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


def boxplot_by_condition(df: pd.DataFrame, axis_col: str, axis_label: str, outpath: Path) -> None:
    order = ["mgmt_non_methylated", "mgmt_methylated"]
    data = [pd.to_numeric(df.loc[df["condition_clean"] == d, axis_col], errors="coerce").dropna() for d in order]

    plt.figure(figsize=(7, 5))
    plt.boxplot(data, tick_labels=order)
    plt.xlabel("Curated MGMT condition")
    plt.ylabel(axis_label)
    plt.title(f"{axis_label} by curated MGMT condition")
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


def barplot_domain_condition(df: pd.DataFrame, outpath: Path) -> None:
    counts = (
        df.groupby(["state_domain", "condition_clean"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=["mgmt_non_methylated", "mgmt_methylated"], fill_value=0)
        .sort_index()
    )

    ax = counts.plot(kind="bar", figsize=(7, 5))
    ax.set_xlabel("State domain")
    ax.set_ylabel("Count")
    ax.set_title("Curated MGMT counts by state domain")
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--state-domains-csv", required=True)
    ap.add_argument("--sample-annotations-curated", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    ensure_outdir(outdir)

    domains = pd.read_csv(args.state_domains_csv)
    ann = load_table(args.sample_annotations_curated)

    merged = domains.merge(ann, on="sample_id", how="left", suffixes=("_dom", "_ann"))

    axis_map = {axis: resolve_axis_column(merged, axis) for axis in AXES}
    cond_col = choose_condition_column(merged)
    merged["condition_clean"] = merged[cond_col].map(normalize_condition_value)

    curated = merged.loc[merged["condition_clean"].isin(["mgmt_methylated", "mgmt_non_methylated"])].copy()

    for axis in AXES:
        axis_col = axis_map[axis]
        boxplot_by_domain(merged, axis_col, axis, outdir / f"{axis.lower()}_by_state_domain.png")
        if len(curated) > 0:
            boxplot_by_condition(curated, axis_col, axis, outdir / f"{axis.lower()}_by_mgmt_condition.png")

    if len(curated) > 0:
        barplot_domain_condition(curated, outdir / "mgmt_counts_by_state_domain.png")

    print(f"[ok] wrote plots to {outdir}")
    print(f"[info] axis column map: {axis_map}")
    print(f"[info] condition column used: {cond_col}")


if __name__ == "__main__":
    main()
