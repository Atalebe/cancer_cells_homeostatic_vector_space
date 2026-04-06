#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def main() -> None:
    cfg = read_config()
    meta_dir = REPO_ROOT / cfg["metadata_dir"]
    out_dir = REPO_ROOT / cfg["results_dir"] / "audit_reports"
    out_dir.mkdir(parents=True, exist_ok=True)

    registry_path = meta_dir / "cell_metadata_registry.parquet"
    if not registry_path.exists():
        raise FileNotFoundError(registry_path)

    reg = pd.read_parquet(registry_path)

    summary = {
        "n_registry_rows": int(len(reg)),
        "n_registry_cols": int(len(reg.columns)),
        "columns": reg.columns.tolist(),
        "n_missing_gsm": int(reg["gsm"].isna().sum()) if "gsm" in reg.columns else None,
        "n_unique_gsm": int(reg["gsm"].nunique()) if "gsm" in reg.columns else None,
        "duplicate_gsm_count": int(reg["gsm"].duplicated().sum()) if "gsm" in reg.columns else None,
    }

    out_json = out_dir / "dataset_registry_audit.json"
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"[ok] wrote {out_json}")


if __name__ == "__main__":
    main()
