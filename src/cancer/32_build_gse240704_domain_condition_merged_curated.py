#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    raise SystemExit(f"unsupported file type: {path}")


def choose_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--state-table", required=True)
    ap.add_argument("--sample-annotations-curated", required=True)
    ap.add_argument("--state-domains-csv", required=True)
    ap.add_argument("--stable-shell-csv", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    state_tbl = read_table(Path(args.state_table)).copy()
    ann = read_table(Path(args.sample_annotations_curated)).copy()
    dom = pd.read_csv(args.state_domains_csv).copy()
    shell = pd.read_csv(args.stable_shell_csv).copy()

    sample_col_state = choose_col(state_tbl, ["sample_id"])
    sample_col_ann = choose_col(ann, ["sample_id"])
    sample_col_dom = choose_col(dom, ["sample_id"])
    sample_col_shell = choose_col(shell, ["sample_id"])

    if not all([sample_col_state, sample_col_ann, sample_col_dom, sample_col_shell]):
        raise SystemExit("sample_id column missing in one or more input tables")

    state_tbl = state_tbl.rename(columns={sample_col_state: "sample_id"})
    ann = ann.rename(columns={sample_col_ann: "sample_id"})
    dom = dom.rename(columns={sample_col_dom: "sample_id"})
    shell = shell.rename(columns={sample_col_shell: "sample_id"})

    merged = state_tbl.merge(
        ann,
        on="sample_id",
        how="left",
        suffixes=("_state_tbl", "_ann"),
    )

    dom_keep = [c for c in ["sample_id", "state_domain"] if c in dom.columns]
    merged = merged.merge(dom[dom_keep], on="sample_id", how="left")

    shell_flag = shell[["sample_id"]].copy()
    shell_flag["candidate_stable_shell"] = True
    merged = merged.merge(shell_flag, on="sample_id", how="left")
    merged["candidate_stable_shell"] = merged["candidate_stable_shell"].fillna(False)

    axis_map = {}
    for axis in ["H", "S", "M", "R", "phi"]:
        candidates = [
            axis,
            f"{axis}_state",
            f"{axis}_state_tbl",
            f"{axis}_ann",
        ]
        chosen = None
        for c in candidates:
            if c in merged.columns:
                chosen = c
                break
        if chosen is None:
            raise SystemExit(f"could not resolve axis column for {axis}")
        axis_map[axis] = chosen
        merged[f"{axis}_dom"] = pd.to_numeric(merged[chosen], errors="coerce")

    cond_col = choose_col(
        merged,
        [
            "placeholder_condition_cur",
            "biological_condition",
            "tumor_status",
            "_condition_clean",
            "placeholder_condition",
        ],
    )
    if cond_col is None:
        merged["placeholder_condition_cur"] = pd.NA
        cond_col = "placeholder_condition_cur"

    out_csv = outdir / "domain_condition_merged_curated.csv"
    merged.to_csv(out_csv, index=False)

    debug = {
        "rows": int(len(merged)),
        "columns": int(len(merged.columns)),
        "condition_column_used": cond_col,
        "axis_column_map": axis_map,
        "n_state_domain_non_null": int(merged["state_domain"].notna().sum()) if "state_domain" in merged.columns else 0,
        "n_stable_shell_true": int(merged["candidate_stable_shell"].fillna(False).astype(bool).sum()),
    }
    with open(outdir / "domain_condition_merged_curated_summary.json", "w", encoding="utf-8") as f:
        json.dump(debug, f, indent=2)

    print(f"[ok] wrote merged curated table: {out_csv}")
    print("[info] summary:")
    print(f"  rows: {debug['rows']}")
    print(f"  columns: {debug['columns']}")
    print(f"  condition_column_used: {debug['condition_column_used']}")
    print(f"  axis_column_map: {debug['axis_column_map']}")
    print(f"  n_state_domain_non_null: {debug['n_state_domain_non_null']}")
    print(f"  n_stable_shell_true: {debug['n_stable_shell_true']}")


if __name__ == "__main__":
    main()
