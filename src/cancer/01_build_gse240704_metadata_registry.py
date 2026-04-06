from __future__ import annotations

from pathlib import Path

from _gse240704_utils import (
    load_config,
    make_sample_registry,
    parse_non_normalized_columns,
    parse_normalized_columns,
    read_header_tsv_gz,
)

def main(config_path: str) -> None:
    cfg = load_config(config_path)

    norm_path = cfg["files"]["normalized_matrix"]
    raw_path = cfg["files"]["non_normalized_matrix"]
    out_path = cfg["outputs"]["sample_registry"]

    norm_cols = read_header_tsv_gz(norm_path)
    raw_cols = read_header_tsv_gz(raw_path)

    norm_samples, _ = parse_normalized_columns(norm_cols)
    raw_samples, _ = parse_non_normalized_columns(raw_cols)

    registry = make_sample_registry(norm_samples, raw_samples)

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    registry.to_csv(out_path, index=False)

    print(f"[ok] wrote sample registry: {out_path}")
    print(f"[info] normalized sample count: {len(norm_samples)}")
    print(f"[info] non-normalized sample count: {len(raw_samples)}")
    print(registry.head(10).to_string(index=False))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    args = parser.parse_args()
    main(args.config)
