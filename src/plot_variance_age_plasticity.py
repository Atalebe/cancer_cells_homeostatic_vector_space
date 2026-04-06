from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REPO = Path.cwd()
RESULTS = REPO / "results" / "gse161895"

def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)

def load_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path)

def save_json(obj: dict, path: Path) -> None:
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)

def rank01(values: pd.Series) -> pd.Series:
    vals = pd.to_numeric(values, errors="coerce")
    vals = vals.fillna(vals.median() if vals.notna().sum() else 0)
    if vals.nunique(dropna=True) <= 1:
        return pd.Series(np.full(len(vals), 0.5), index=vals.index)
    return vals.rank(method="average", pct=True)

def style_axes(ax, title: str, xlabel: str = "", ylabel: str = "") -> None:
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25)

def main() -> None:
    outdir = RESULTS / "visual_package" / "variance_age_plasticity"
    ensure_dir(outdir / "plots")
    vs = load_csv(RESULTS / "d2_1_variance_scaling" / "variance_scaling_summary.csv")
    age = load_csv(RESULTS / "age_axis_proxy" / "state_age_proxy_tests.csv")
    plast = load_csv(RESULTS / "plasticity_reversion_diagnostics" / "plasticity_reversion_domain_tests.csv")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(vs["group"], vs["slope_logvar_vs_logmean"])
    style_axes(ax, "Variance scaling slopes across key compartments", "", "Slope log(var) vs log(mean)")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(outdir / "plots" / "variance_scaling_slopes.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    labels = [f"{a} vs {b}" for a, b in zip(age["group_a"], age["group_b"])]
    ax.bar(labels, age["delta_median_a_minus_b"])
    style_axes(ax, "State-age proxy contrasts", "", "Delta median")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(outdir / "plots" / "state_age_proxy_contrasts.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(plast["metric"], plast["delta_median_a_minus_b"])
    style_axes(ax, "Plasticity and donor reversion contrasts, D1 vs D2", "", "Delta median D1 minus D2")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(outdir / "plots" / "plasticity_reversion_contrasts.png", dpi=220)
    plt.close(fig)

if __name__ == "__main__":
    main()
