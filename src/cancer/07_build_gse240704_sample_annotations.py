from __future__ import annotations

from pathlib import Path

import pandas as pd

from _gse240704_utils import load_config


def main(config_path: str) -> None:
    cfg = load_config(config_path)

    registry_path = "data/metadata/gse240704/sample_registry.csv"
    proxy_path = cfg["outputs"]["hrsm_proxy_table"]
    state_path = cfg["outputs"]["state_table"]
    out_path = "data/processed/gse240704/sample_annotations.parquet"
    out_csv = "data/processed/gse240704/sample_annotations.csv"

    registry = pd.read_csv(registry_path)
    proxy = pd.read_csv(proxy_path)
    state = pd.read_parquet(state_path)

    df = registry.merge(proxy, on="sample_id", how="left", suffixes=("", "_proxy"))
    df = df.merge(state, on="sample_id", how="left", suffixes=("", "_state"))

    if "phi" in df.columns:
        df["phi_rank_desc"] = df["phi"].rank(method="dense", ascending=False)

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, index=False)
    df.to_csv(out_csv, index=False)

    print(f"[ok] wrote merged sample annotations: {out_path}")
    print(f"[info] rows: {len(df)}")
    print(f"[info] columns: {len(df.columns)}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    args = parser.parse_args()
    main(args.config)
