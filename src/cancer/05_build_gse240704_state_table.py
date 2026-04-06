from __future__ import annotations

from pathlib import Path

import pandas as pd

from _gse240704_utils import load_config

def main(config_path: str) -> None:
    cfg = load_config(config_path)

    registry_path = cfg["outputs"]["sample_registry"]
    proxy_path = cfg["outputs"]["hrsm_proxy_table"]
    out_path = cfg["outputs"]["state_table"]

    registry = pd.read_csv(registry_path)
    proxy = pd.read_csv(proxy_path)

    state = registry.merge(proxy, on="sample_id", how="left")

    state["state_available"] = state["H"].notna() & state["S"].notna() & state["M"].notna() & state["R"].notna()

    cols = [
        "sample_id",
        "sample_number",
        "in_normalized_matrix",
        "in_non_normalized_matrix",
        "placeholder_condition",
        "placeholder_group",
        "placeholder_batch",
        "placeholder_note",
        "state_available",
        "H",
        "R",
        "S",
        "M",
        "phi",
        "H_raw",
        "R_raw",
        "S_raw",
        "M_raw",
        "mean_neighbor_distance",
        "pca_explained_variance_pc1",
        "pca_explained_variance_pc2",
    ]
    state = state[cols]

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    state.to_parquet(out_path, index=False)

    print(f"[ok] wrote state table: {out_path}")
    print(f"[info] rows: {len(state)}")
    print(f"[info] state_available: {int(state['state_available'].sum())}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    args = parser.parse_args()
    main(args.config)
