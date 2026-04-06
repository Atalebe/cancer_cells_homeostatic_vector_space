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
    outdir = RESULTS / "visual_package" / "trajectory_visuals"
    ensure_dir(outdir / "plots")
    cent = load_csv(RESULTS / "transition_trajectory_pass" / "trajectory_centroids.csv")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(cent["H"], cent["R"], marker="o")
    for _, row in cent.iterrows():
        ax.annotate(row["trajectory_group"], (row["H"], row["R"]), fontsize=9)
    style_axes(ax, "Formal trajectory pass, H-R plane", "H", "R")
    plt.tight_layout()
    plt.savefig(outdir / "plots" / "trajectory_hr_plane.png", dpi=220)
    plt.close(fig)

    metrics = ["H", "S", "M", "R"]
    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(metrics))
    for _, row in cent.iterrows():
        ax.plot(x, [row[m] for m in metrics], marker="o", label=row["trajectory_group"])
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.legend()
    style_axes(ax, "Trajectory centroid profiles across HRSM", "", "Centroid value")
    plt.tight_layout()
    plt.savefig(outdir / "plots" / "trajectory_centroid_profiles.png", dpi=220)
    plt.close(fig)

if __name__ == "__main__":
    main()
