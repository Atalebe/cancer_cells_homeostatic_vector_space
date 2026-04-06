from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


AXES = ["H", "S", "M", "R", "phi"]
VALID_CONDITIONS = {"mgmt_methylated", "mgmt_non_methylated"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run curated MGMT condition-axis tests and stable-shell overlay for GSE240704."
    )
    parser.add_argument("--state-table", required=True)
    parser.add_argument("--sample-annotations-curated", required=True)
    parser.add_argument("--stable-shell-csv", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument(
        "--state-domain-table",
        default="results/gse240704/state_domains/state_domain_assignments.csv",
        help="Optional CSV with sample_id,state_domain columns.",
    )
    return parser.parse_args()


def load_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Missing required input: {path}")
    if p.suffix.lower() == ".parquet":
        return pd.read_parquet(p)
    return pd.read_csv(p)


def effect_size_rank_biserial(x: pd.Series, y: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce").dropna().to_numpy()
    y = pd.to_numeric(y, errors="coerce").dropna().to_numpy()
    if len(x) == 0 or len(y) == 0:
        return np.nan
    u = mannwhitneyu(x, y, alternative="two-sided").statistic
    return 2.0 * u / (len(x) * len(y)) - 1.0


def maybe_attach_state_domain(df: pd.DataFrame, state_domain_table: str) -> pd.DataFrame:
    if "state_domain" in df.columns:
        return df

    p = Path(state_domain_table)
    if not p.exists():
        df["state_domain"] = pd.NA
        return df

    dom = pd.read_csv(p)
    if not {"sample_id", "state_domain"}.issubset(dom.columns):
        df["state_domain"] = pd.NA
        return df

    dom = dom[["sample_id", "state_domain"]].drop_duplicates(subset=["sample_id"])
    return df.merge(dom, on="sample_id", how="left")


def choose_condition_column(df: pd.DataFrame) -> str | None:
    candidates = [
        "placeholder_condition_cur",
        "biological_condition",
        "placeholder_condition",
        "placeholder_condition_ann",
        "tumor_status",
    ]
    best_col = None
    best_hits = -1

    for col in candidates:
        if col not in df.columns:
            continue
        vals = set(
            df[col]
            .dropna()
            .astype(str)
            .str.strip()
            .str.lower()
            .tolist()
        )
        hits = len(vals & VALID_CONDITIONS)
        if hits > best_hits:
            best_hits = hits
            best_col = col

    return best_col


def ensure_axis_columns(df: pd.DataFrame) -> pd.DataFrame:
    for axis in AXES:
        if axis in df.columns:
            continue

        candidates = [
            f"{axis}_state",
            f"{axis}_x",
            f"{axis}_ann",
            f"{axis}_y",
        ]
        found = next((c for c in candidates if c in df.columns), None)
        if found is not None:
            df[axis] = df[found]

    return df


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    st = load_table(args.state_table)
    ann = load_table(args.sample_annotations_curated)

    df = st.merge(ann, on="sample_id", how="left", suffixes=("_state_tbl", "_ann"))
    df = maybe_attach_state_domain(df, args.state_domain_table)
    df = ensure_axis_columns(df)

    condition_col = choose_condition_column(df)
    if condition_col is None:
        raise ValueError("No usable condition column found in merged dataframe.")

    df["_condition_clean"] = (
        df[condition_col]
        .astype("string")
        .str.strip()
        .str.lower()
    )

    sub = df[df["_condition_clean"].isin(VALID_CONDITIONS)].copy()

    debug = {
        "merged_shape": list(df.shape),
        "condition_column_used": condition_col,
        "available_condition_like_columns": [
            c for c in df.columns if "condition" in c.lower() or "tumor" in c.lower()
        ],
        "condition_value_counts_top": df[condition_col]
        .astype("string")
        .fillna("<NA>")
        .value_counts(dropna=False)
        .head(20)
        .to_dict(),
        "subset_shape": list(sub.shape),
        "axis_columns_present": {axis: (axis in sub.columns) for axis in AXES},
    }

    with open(outdir / "condition_axis_debug.json", "w", encoding="utf-8") as f:
        json.dump(debug, f, indent=2)

    rows = []
    a = sub[sub["_condition_clean"] == "mgmt_methylated"]
    b = sub[sub["_condition_clean"] == "mgmt_non_methylated"]

    for axis in AXES:
        if axis not in sub.columns:
            continue

        xa = pd.to_numeric(a[axis], errors="coerce").dropna()
        xb = pd.to_numeric(b[axis], errors="coerce").dropna()

        if len(xa) == 0 or len(xb) == 0:
            continue

        stat = mannwhitneyu(xa, xb, alternative="two-sided")
        rows.append(
            {
                "condition_a": "mgmt_methylated",
                "condition_b": "mgmt_non_methylated",
                "axis": axis,
                "n_a": int(len(xa)),
                "n_b": int(len(xb)),
                "median_a": float(np.nanmedian(xa)),
                "median_b": float(np.nanmedian(xb)),
                "delta_median_a_minus_b": float(np.nanmedian(xa) - np.nanmedian(xb)),
                "u_statistic": float(stat.statistic),
                "p_value": float(stat.pvalue),
                "rank_biserial": float(effect_size_rank_biserial(xa, xb)),
            }
        )

    if rows:
        out_tests = pd.DataFrame(rows).sort_values("axis").reset_index(drop=True)
    else:
        out_tests = pd.DataFrame(
            columns=[
                "condition_a",
                "condition_b",
                "axis",
                "n_a",
                "n_b",
                "median_a",
                "median_b",
                "delta_median_a_minus_b",
                "u_statistic",
                "p_value",
                "rank_biserial",
            ]
        )

    out_tests.to_csv(outdir / "condition_axis_tests_curated.csv", index=False)

    stable = pd.read_csv(args.stable_shell_csv)
    stable = stable[["sample_id"]].drop_duplicates().copy()

    overlay_cols = [
        "sample_id",
        condition_col,
        "tumor_status",
        "characteristics_ch1",
        "H",
        "S",
        "M",
        "R",
        "phi",
        "state_domain",
    ]
    overlay_cols = [c for c in overlay_cols if c in df.columns]

    stable_overlay = stable.merge(
        df[overlay_cols].drop_duplicates(subset=["sample_id"]),
        on="sample_id",
        how="left",
    )
    stable_overlay.to_csv(outdir / "stable_shell_mgmt_overlay.csv", index=False)

    summary = {
        "n_total_merged": int(len(df)),
        "condition_column_used": condition_col,
        "n_total_curated_subset": int(len(sub)),
        "n_mgmt_methylated": int((sub["_condition_clean"] == "mgmt_methylated").sum()),
        "n_mgmt_non_methylated": int((sub["_condition_clean"] == "mgmt_non_methylated").sum()),
        "n_stable_shell_samples": int(len(stable_overlay)),
        "n_axes_tested": int(len(out_tests)),
    }

    with open(outdir / "condition_axis_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote curated condition-axis tests:", outdir / "condition_axis_tests_curated.csv")
    print("[ok] wrote stable-shell MGMT overlay:", outdir / "stable_shell_mgmt_overlay.csv")
    print("[ok] wrote debug info:", outdir / "condition_axis_debug.json")
    print("[ok] wrote summary:", outdir / "condition_axis_summary.json")
    print("[info] summary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
