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
    outdir = RESULTS / "visual_package" / "program_shift_panels"
    ensure_dir(outdir / "plots")

    d2 = load_csv(RESULTS / "d2_1_treatment_response" / "d2_1_treatment_program_overlay_tests.csv")
    p3 = load_csv(RESULTS / "patient3_d2_1_treatment_sensitivity" / "patient3_d2_1_treatment_program_overlay_tests.csv")
    metrics = ["mito_fraction", "mitochondrial_encoded_score", "immune_leaning_score", "ribosomal_translation_score", "antigen_presentation_score"]

    d2 = d2[d2["metric"].isin(metrics)].copy().set_index("metric")
    p3 = p3[p3["metric"].isin(metrics)].copy().set_index("metric")

    x = np.arange(len(metrics))
    width = 0.38
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(x - width/2, d2.reindex(metrics)["delta_median_a_minus_b"], width=width, label="D2_1 treated vs untreated")
    ax.bar(x + width/2, p3.reindex(metrics)["delta_median_a_minus_b"], width=width, label="Patient3 D2_1 treated vs untreated")
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, rotation=30, ha="right")
    ax.legend()
    style_axes(ax, "Treatment program shifts, whole D2_1 and Patient3 sensitivity", "", "Delta median treated minus untreated")
    plt.tight_layout()
    plt.savefig(outdir / "plots" / "treatment_program_shift_comparison.png", dpi=220)
    plt.close(fig)

if __name__ == "__main__":
    main()
