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
    outdir = RESULTS / "visual_package" / "go_mechanism_panels"
    ensure_dir(outdir / "plots")
    up = load_csv(RESULTS / "d2_1_r_go_enrichment" / "up_go_enrichment_results.csv")
    down = load_csv(RESULTS / "d2_1_r_go_enrichment" / "down_go_enrichment_results.csv")

    for df in (up, down):
        df["Adjusted P-value"] = pd.to_numeric(df["Adjusted P-value"], errors="coerce").clip(lower=1e-300)
        df["neglog10_adj_p"] = -np.log10(df["Adjusted P-value"])

    for name, df in [("up", up), ("down", down)]:
        top = df.nsmallest(12, "Adjusted P-value").sort_values("neglog10_adj_p")
        fig, ax = plt.subplots(figsize=(10, 7))
        ax.barh(range(len(top)), top["neglog10_adj_p"].values)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top["Term"])
        style_axes(ax, f"Top GO terms, {name.replace('up','treated-up/R-positive').replace('down','treated-down/R-negative')}", "-log10 adjusted p", "")
        plt.tight_layout()
        plt.savefig(outdir / "plots" / f"go_{name}_top_terms.png", dpi=220)
        plt.close(fig)

if __name__ == "__main__":
    main()
