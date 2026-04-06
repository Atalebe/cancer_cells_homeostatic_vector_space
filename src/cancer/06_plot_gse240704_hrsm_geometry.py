from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from _gse240704_utils import load_config


def _scatter(df: pd.DataFrame, x: str, y: str, out_path: str) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    ax.scatter(df[x], df[y], s=22, alpha=0.75)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_title(f"{x} vs {y}")
    ax.axhline(0.0, linewidth=0.8)
    ax.axvline(0.0, linewidth=0.8)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def main(config_path: str) -> None:
    cfg = load_config(config_path)

    state_path = cfg["outputs"]["state_table"]
    out_dir = Path("results/gse240704/geometry")
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(state_path)

    needed = ["sample_id", "H", "S", "M", "R", "phi"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise RuntimeError(f"Missing required columns in state table: {missing}")

    _scatter(df, "H", "S", str(out_dir / "hs_plane.png"))
    _scatter(df, "H", "M", str(out_dir / "hm_plane.png"))
    _scatter(df, "S", "R", str(out_dir / "sr_plane.png"))
    _scatter(df, "M", "R", str(out_dir / "mr_plane.png"))
    _scatter(df, "phi", "R", str(out_dir / "phi_r_plane.png"))

    summary = {
        "n_samples": int(len(df)),
        "phi_min": float(df["phi"].min()),
        "phi_median": float(df["phi"].median()),
        "phi_max": float(df["phi"].max()),
        "H_median": float(df["H"].median()),
        "S_median": float(df["S"].median()),
        "M_median": float(df["M"].median()),
        "R_median": float(df["R"].median()),
    }
    pd.Series(summary).to_csv(out_dir / "geometry_summary.csv", header=["value"])

    print(f"[ok] wrote geometry plots to {out_dir}")
    print("[info] summary")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    args = parser.parse_args()
    main(args.config)
