from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from _gse240704_utils import (
    count_lines_gz,
    load_config,
    parse_non_normalized_columns,
    parse_normalized_columns,
    read_header_tsv_gz,
)

def audit_normalized(cfg: dict) -> dict:
    path = cfg["files"]["normalized_matrix"]
    det_thr = cfg["qc"]["detection_pval_threshold"]

    cols = read_header_tsv_gz(path)
    sample_ids, _ = parse_normalized_columns(cols)
    n_rows = count_lines_gz(path) - 1

    fail_counts = {s: 0 for s in sample_ids}
    miss_counts = {s: 0 for s in sample_ids}
    total_rows = 0

    usecols = ["ID_REF"]
    for s in sample_ids:
        usecols.extend([s, "Detection Pval"])

    chunk_iter = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        chunksize=5000,
        low_memory=False,
    )

    for chunk in chunk_iter:
        total_rows += len(chunk)
        pairs = chunk.columns.tolist()[1:]
        pair_idx = 0
        for s in sample_ids:
            value_col = pairs[pair_idx]
            det_col = pairs[pair_idx + 1]
            pair_idx += 2
            vals = pd.to_numeric(chunk[value_col], errors="coerce")
            dets = pd.to_numeric(chunk[det_col], errors="coerce")
            fail_counts[s] += int((dets > det_thr).fillna(False).sum())
            miss_counts[s] += int(vals.isna().sum() + dets.isna().sum())

    return {
        "file": path,
        "n_probes_estimated": int(n_rows),
        "n_samples": len(sample_ids),
        "sample_ids_head": sample_ids[:10],
        "detection_pval_threshold": det_thr,
        "per_sample_bad_probe_fraction": {
            s: fail_counts[s] / max(total_rows, 1) for s in sample_ids
        },
        "per_sample_missing_entry_fraction": {
            s: miss_counts[s] / max(total_rows * 2, 1) for s in sample_ids
        },
    }


def audit_non_normalized(cfg: dict) -> dict:
    path = cfg["files"]["non_normalized_matrix"]
    cols = read_header_tsv_gz(path)
    sample_ids, mapping = parse_non_normalized_columns(cols)
    n_rows = count_lines_gz(path) - 1

    completeness = {}
    for s in sample_ids:
        completeness[s] = {
            "has_u_col": "u_col" in mapping[s],
            "has_m_col": "m_col" in mapping[s],
            "has_det_col": "det_col" in mapping[s],
        }

    return {
        "file": path,
        "n_probes_estimated": int(n_rows),
        "n_samples": len(sample_ids),
        "sample_ids_head": sample_ids[:10],
        "triplet_completeness": completeness,
    }


def main(config_path: str) -> None:
    cfg = load_config(config_path)
    out_path = cfg["outputs"]["matrix_audit"]

    normalized = audit_normalized(cfg)
    non_normalized = audit_non_normalized(cfg)

    bad_probe_warn = cfg["qc"]["bad_probe_fraction_warn"]
    missing_warn = cfg["qc"]["missing_fraction_warn"]

    flagged_samples = []
    for s, frac in normalized["per_sample_bad_probe_fraction"].items():
        miss_frac = normalized["per_sample_missing_entry_fraction"][s]
        if frac > bad_probe_warn or miss_frac > missing_warn:
            flagged_samples.append(
                {
                    "sample_id": s,
                    "bad_probe_fraction": frac,
                    "missing_entry_fraction": miss_frac,
                }
            )

    payload = {
        "dataset_id": cfg["dataset_id"],
        "normalized": normalized,
        "non_normalized": non_normalized,
        "warn_thresholds": {
            "bad_probe_fraction_warn": bad_probe_warn,
            "missing_fraction_warn": missing_warn,
        },
        "flagged_samples": flagged_samples,
    }

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    print(f"[ok] wrote matrix audit: {out_path}")
    print(f"[info] normalized samples: {normalized['n_samples']}")
    print(f"[info] non-normalized samples: {non_normalized['n_samples']}")
    print(f"[info] flagged samples: {len(flagged_samples)}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    args = parser.parse_args()
    main(args.config)
