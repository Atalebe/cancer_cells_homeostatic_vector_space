from __future__ import annotations

from pathlib import Path

import pandas as pd

from _gse240704_utils import load_config, parse_normalized_columns, read_header_tsv_gz, safe_float_frame


def main(config_path: str) -> None:
    cfg = load_config(config_path)

    in_path = cfg["files"]["normalized_matrix"]
    beta_out = cfg["outputs"]["beta_matrix"]
    det_out = cfg["outputs"]["detection_pval_matrix"]
    det_thr = cfg["qc"]["detection_pval_threshold"]

    cols = read_header_tsv_gz(in_path)
    sample_ids, _ = parse_normalized_columns(cols)

    beta_chunks = []
    det_chunks = []

    reader = pd.read_csv(
        in_path,
        sep="\t",
        compression="gzip",
        chunksize=5000,
        low_memory=False,
    )

    for chunk_idx, chunk in enumerate(reader, start=1):
        pair_cols = chunk.columns.tolist()[1:]

        value_cols = []
        det_cols = []
        pair_idx = 0
        for _ in sample_ids:
            value_cols.append(pair_cols[pair_idx])
            det_cols.append(pair_cols[pair_idx + 1])
            pair_idx += 2

        beta = chunk[value_cols].copy()
        beta.columns = sample_ids
        beta = safe_float_frame(beta)

        det = chunk[det_cols].copy()
        det.columns = sample_ids
        det = safe_float_frame(det)

        beta = beta.mask(det > det_thr)

        beta = pd.concat([chunk[["ID_REF"]].copy(), beta], axis=1)
        det = pd.concat([chunk[["ID_REF"]].copy(), det], axis=1)

        beta_chunks.append(beta)
        det_chunks.append(det)

        print(f"[info] processed chunk {chunk_idx}")

    beta_df = pd.concat(beta_chunks, axis=0, ignore_index=True)
    det_df = pd.concat(det_chunks, axis=0, ignore_index=True)

    Path(beta_out).parent.mkdir(parents=True, exist_ok=True)
    beta_df.to_parquet(beta_out, index=False)
    det_df.to_parquet(det_out, index=False)

    print(f"[ok] wrote beta matrix: {beta_out}")
    print(f"[ok] wrote detection p-value matrix: {det_out}")
    print(f"[info] beta matrix shape: {beta_df.shape}")
    print(f"[info] detection matrix shape: {det_df.shape}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    args = parser.parse_args()
    main(args.config)
