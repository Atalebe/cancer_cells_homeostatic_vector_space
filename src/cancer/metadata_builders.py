from __future__ import annotations

from pathlib import Path
import yaml
import pandas as pd

from common.io import read_table, write_table


def build_metadata_registry(config_path: str | Path, output_path: str | Path) -> None:
    with open(config_path, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    metadata_file = cfg["paths"]["metadata_file"]
    barcode_col = cfg["metadata"]["barcode_column"]
    metadata = read_table(metadata_file)
    metadata[barcode_col] = metadata[barcode_col].astype(str)
    metadata = metadata.drop_duplicates(subset=[barcode_col]).copy()
    write_table(metadata, output_path)
