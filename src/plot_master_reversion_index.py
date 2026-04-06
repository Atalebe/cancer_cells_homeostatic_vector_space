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
    outdir = RESULTS / "visual_package" / "reversion_index_master"
    ensure_dir(outdir / "plots")

    bridge = load_csv(RESULTS / "d2_1_r_treatment_bridge" / "d2_1_r_treatment_bridge_up_hits.csv")
    rev = load_csv(RESULTS / "reversion_candidate_scan" / "reversion_candidate_full_scan.csv")

    for df in (bridge, rev):
        if "gene_id_stripped" not in df.columns and "gene_id" in df.columns:
            df["gene_id_stripped"] = df["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
        if "gene_symbol" not in df.columns:
            df["gene_symbol"] = ""

    merge_cols = [c for c in ["gene_id_stripped", "gene_symbol", "bridge_score", "bridge_direction"] if c in bridge.columns]
    merged = rev.merge(bridge[merge_cols].drop_duplicates("gene_id_stripped"), on="gene_id_stripped", how="left", suffixes=("", "_bridge"))
    if "gene_symbol_bridge" in merged.columns:
        merged["gene_symbol"] = merged["gene_symbol"].replace("", np.nan).fillna(merged["gene_symbol_bridge"]).fillna("")

    candidate_cols = [c for c in [
        "toward_donor_score",
        "reversion_score",
        "delta_treated_minus_untreated",
        "delta_log1p_treated_minus_untreated",
        "bridge_score",
        "donor_reversion_score",
    ] if c in merged.columns]

    if not candidate_cols:
        # fall back to numeric score-like columns
        candidate_cols = [c for c in merged.columns if any(k in c.lower() for k in ["score", "delta", "reversion", "donor"])]

    composite = pd.Series(0.0, index=merged.index)
    used = 0
    for c in candidate_cols[:4]:
        s = pd.to_numeric(merged[c], errors="coerce")
        if s.notna().sum():
            composite += rank01(s)
            used += 1
    if used == 0:
        raise ValueError("No usable numeric score columns found for master reversion index.")
    merged["master_reversion_index"] = composite / used

    if "bridge_score" in merged.columns:
        merged["biophysical_vector_index"] = 0.6 * merged["master_reversion_index"] + 0.4 * rank01(merged["bridge_score"])
    else:
        merged["biophysical_vector_index"] = merged["master_reversion_index"]

    cols = [c for c in ["gene_id_stripped", "gene_symbol", "master_reversion_index", "biophysical_vector_index", "bridge_score", "bridge_direction"] + candidate_cols[:4] if c in merged.columns]
    summary = merged[cols].sort_values("biophysical_vector_index", ascending=False).copy()
    summary.to_csv(outdir / "reversion_index_master_table.csv", index=False)

    top = summary.head(25).copy()
    labels = top["gene_symbol"].replace("", np.nan).fillna(top["gene_id_stripped"])

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.barh(range(len(top)), top["biophysical_vector_index"].values)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    style_axes(ax, "Master reversion index, donor-directed biophysical response genes", "Composite index", "")
    plt.tight_layout()
    plt.savefig(outdir / "plots" / "master_reversion_index_top25.png", dpi=220)
    plt.close(fig)

    if "bridge_score" in summary.columns:
        fig, ax = plt.subplots(figsize=(8, 7))
        ax.scatter(summary["master_reversion_index"], summary["bridge_score"], alpha=0.5)
        for _, row in summary.head(15).iterrows():
            label = row["gene_symbol"] if str(row.get("gene_symbol", "")) else row["gene_id_stripped"]
            ax.annotate(label, (row["master_reversion_index"], row["bridge_score"]), fontsize=8)
        style_axes(ax, "Biophysical vector of treatment response", "Master reversion index", "R-bridge score")
        plt.tight_layout()
        plt.savefig(outdir / "plots" / "reversion_vs_bridge_scatter.png", dpi=220)
        plt.close(fig)

    save_json({
        "n_genes": int(len(summary)),
        "top_labels": labels.head(10).tolist(),
        "score_columns_used": candidate_cols[:4],
    }, outdir / "reversion_index_master_summary.json")

if __name__ == "__main__":
    main()
