#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def main() -> None:
    cfg = read_config()
    proc_dir = REPO_ROOT / cfg["processed_dir"]

    proxy_df = pd.read_csv(proc_dir / "hrsm_proxy_table.csv")
    proxy_df["phi"] = np.sqrt((proxy_df[["H", "S", "M", "R"]] ** 2).sum(axis=1))
    proxy_df.to_parquet(proc_dir / "state_table.parquet", index=False)
    proxy_df.to_csv(proc_dir / "state_table.csv", index=False)

    print("[ok] wrote state table")


if __name__ == "__main__":
    main()
