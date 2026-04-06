from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy.stats import mannwhitneyu

from _gse240704_utils import load_config


def _safe_test(df: pd.DataFrame, group_col: str, value_col: str):
    sub = df[[group_col, value_col]].dropna().copy()
    groups = list(sub[group_col].astype(str).unique())
    if len(groups) != 2:
        return None
    a = sub.loc[sub[group_col].astype(str) == groups[0], value_col]
    b = sub.loc[sub[group_col].astype(str) == groups[1], value_col]
    if len(a) < 3 or len(b) < 3:
        return None
    stat, p = mannwhitneyu(a, b, alternative="two-sided")
    return {
        "group_col": group_col,
        "value_col": value_col,
        "group_a": groups[0],
        "group_b": groups[1],
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(a.median()),
        "median_b": float(b.median()),
        "u_stat": float(stat),
        "p_value": float(p),
    }


def main(config_path: str) -> None:
    cfg = load_config(config_path)

    ann_path = "data/processed/gse240704/sample_annotations.parquet"
    out_dir = Path("results/gse240704/condition_axis_tests")
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(ann_path)

    candidate_group_cols = [
        "placeholder_condition",
        "placeholder_group",
        "placeholder_batch",
    ]
    value_cols = ["H", "S", "M", "R", "phi"]

    rows = []
    for gcol in candidate_group_cols:
        if gcol not in df.columns:
            continue
        non_null = df[gcol].dropna()
        if non_null.empty:
            print(f"[info] skipping {gcol}, no annotation yet")
            continue
        for vcol in value_cols:
            if vcol not in df.columns:
                continue
            result = _safe_test(df, gcol, vcol)
            if result is not None:
                rows.append(result)

    out = pd.DataFrame(rows)
    out.to_csv(out_dir / "condition_axis_tests.csv", index=False)

    print(f"[ok] wrote condition axis test table: {out_dir / 'condition_axis_tests.csv'}")
    print(f"[info] rows: {len(out)}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    args = parser.parse_args()
    main(args.config)
