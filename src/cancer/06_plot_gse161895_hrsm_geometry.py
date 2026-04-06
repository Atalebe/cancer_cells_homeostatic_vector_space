#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import yaml


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def save_scatter(df: pd.DataFrame, x: str, y: str, outpath: Path) -> None:
    plt.figure(figsize=(6, 5))
    plt.scatter(df[x], df[y], s=10, alpha=0.7)
    plt.xlabel(x)
    plt.ylabel(y)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def main() -> None:
    cfg = read_config()
    proc_dir = REPO_ROOT / cfg["processed_dir"]
    fig_dir = REPO_ROOT / cfg["results_dir"] / "geometry"
    fig_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(proc_dir / "state_table.parquet")
    save_scatter(df, "H", "S", fig_dir / "hs_plane.png")
    save_scatter(df, "H", "M", fig_dir / "hm_plane.png")
    save_scatter(df, "S", "R", fig_dir / "sr_plane.png")
    save_scatter(df, "M", "R", fig_dir / "mr_plane.png")
    save_scatter(df, "phi", "R", fig_dir / "phi_r_plane.png")

    print("[ok] wrote geometry figures")


if __name__ == "__main__":
    main()
