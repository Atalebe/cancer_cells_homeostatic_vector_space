from __future__ import annotations

import gzip
import json
import math
import re
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import yaml


def load_config(path: str | Path) -> Dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def ensure_parent(path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def robust_z(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    med = np.nanmedian(x.values)
    mad = np.nanmedian(np.abs(x.values - med))
    if not np.isfinite(mad) or mad == 0:
        sd = np.nanstd(x.values)
        if not np.isfinite(sd) or sd == 0:
            return pd.Series(np.zeros(len(x)), index=x.index)
        return pd.Series((x.values - med) / sd, index=x.index)
    return pd.Series(0.67448975 * (x.values - med) / mad, index=x.index)


def read_header_tsv_gz(path: str | Path) -> List[str]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
        header = f.readline().rstrip("\n")
    return header.split("\t")


def read_first_n_rows_tsv_gz(path: str | Path, nrows: int = 5) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="gzip", nrows=nrows, dtype=str, low_memory=False)


def count_lines_gz(path: str | Path) -> int:
    n = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
        for _ in f:
            n += 1
    return n


def parse_normalized_columns(columns: List[str]) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    if not columns or columns[0] != "ID_REF":
        raise ValueError("Expected first column to be ID_REF")
    sample_ids: List[str] = []
    mapping: Dict[str, Dict[str, str]] = {}
    i = 1
    while i < len(columns):
        value_col = columns[i]
        if i + 1 >= len(columns):
            raise ValueError(f"Unpaired normalized column at position {i}: {value_col}")
        det_col = columns[i + 1]
        if det_col != "Detection Pval":
            raise ValueError(f"Expected 'Detection Pval' at position {i+1}, found {det_col}")
        sample_id = value_col.strip()
        sample_ids.append(sample_id)
        mapping[sample_id] = {"value_col": value_col, "det_col": det_col}
        i += 2
    return sample_ids, mapping


_TRIPLET_RE = re.compile(r"^(SAMPLE\s+\d+)\s+(Unmethylated Signal|Methylated Signal|Detection Pval)$", flags=re.IGNORECASE)


def parse_non_normalized_columns(columns: List[str]) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    if not columns or columns[0] != "ID_REF":
        raise ValueError("Expected first column to be ID_REF")
    sample_ids: List[str] = []
    mapping: Dict[str, Dict[str, str]] = {}
    for col in columns[1:]:
        m = _TRIPLET_RE.match(col.strip())
        if not m:
            raise ValueError(f"Unexpected non-normalized column format: {col}")
        sample_id = m.group(1).strip()
        kind = m.group(2).strip().lower()
        if sample_id not in mapping:
            mapping[sample_id] = {}
            sample_ids.append(sample_id)
        if "unmethylated" in kind:
            mapping[sample_id]["u_col"] = col
        elif "methylated" in kind:
            mapping[sample_id]["m_col"] = col
        elif "detection" in kind:
            mapping[sample_id]["det_col"] = col
    return sample_ids, mapping


def write_json(obj: Dict, path: str | Path) -> None:
    ensure_parent(path)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)


def safe_float_frame(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out


def beta_from_raw(m: pd.DataFrame, u: pd.DataFrame, offset: float = 100.0) -> pd.DataFrame:
    return m / (m + u + offset)


def m_value_from_beta(beta: pd.DataFrame, eps: float = 1e-6) -> pd.DataFrame:
    clipped = beta.clip(lower=eps, upper=1.0 - eps)
    return np.log2(clipped / (1.0 - clipped))


def entropy_from_distances(distances: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = 1.0 / np.clip(distances, 1e-12, None)
    inv_sum = inv.sum(axis=1, keepdims=True)
    probs = inv / np.clip(inv_sum, 1e-12, None)
    with np.errstate(divide="ignore", invalid="ignore"):
        ent = -(probs * np.log(np.clip(probs, 1e-12, None))).sum(axis=1)
    return ent


def make_sample_registry(
    normalized_samples: List[str],
    non_normalized_samples: List[str],
) -> pd.DataFrame:
    union_samples = sorted(set(normalized_samples) | set(non_normalized_samples), key=lambda s: int(re.findall(r"\d+", s)[0]))
    rows = []
    for s in union_samples:
        sample_num = int(re.findall(r"\d+", s)[0])
        rows.append(
            {
                "sample_id": s,
                "sample_number": sample_num,
                "in_normalized_matrix": s in normalized_samples,
                "in_non_normalized_matrix": s in non_normalized_samples,
                "placeholder_condition": pd.NA,
                "placeholder_group": pd.NA,
                "placeholder_batch": pd.NA,
                "placeholder_note": pd.NA,
            }
        )
    return pd.DataFrame(rows).sort_values("sample_number").reset_index(drop=True)
