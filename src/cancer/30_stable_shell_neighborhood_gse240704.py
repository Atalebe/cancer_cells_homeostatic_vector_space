#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


AXES = ["H", "S", "M", "R"]


def resolve_axis_cols(df: pd.DataFrame) -> dict[str, str]:
    out = {}
    for axis in AXES + ["phi"]:
        candidates = [axis, f"{axis}_dom", f"{axis}_state"]
        chosen = next((c for c in candidates if c in df.columns), None)
        if chosen is None:
            raise KeyError(f"could not resolve source column for axis '{axis}'")
        out[axis] = chosen
    return out


def euclidean(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.sqrt(((a - b) ** 2).sum(axis=1))


def median_or_nan(x: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce").dropna()
    if len(x) == 0:
        return float("nan")
    return float(np.median(x))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--merged-curated-csv", required=True)
    ap.add_argument("--stable-shell-csv", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--k", type=int, default=15)
    args = ap.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.merged_curated_csv)
    shell = pd.read_csv(args.stable_shell_csv)

    if "sample_id" not in df.columns or "sample_id" not in shell.columns:
        raise SystemExit("sample_id column required in both inputs")

    axis_map = resolve_axis_cols(df)

    work = df.copy()
    for axis in AXES:
        work[axis] = pd.to_numeric(work[axis_map[axis]], errors="coerce")

    work = work.dropna(subset=AXES).reset_index(drop=True)

    shell_ids = shell["sample_id"].astype(str).tolist()
    shell_work = work.loc[work["sample_id"].astype(str).isin(shell_ids)].copy()

    if shell_work.empty:
        raise SystemExit("no stable shell samples found in merged curated table")

    neighbor_rows = []
    summary_rows = []

    space = work[AXES].to_numpy(dtype=float)

    for _, row in shell_work.iterrows():
        sid = str(row["sample_id"])
        anchor = row[AXES].to_numpy(dtype=float).reshape(1, -1)
        dist = euclidean(space, anchor)
        tmp = work.copy()
        tmp["distance_to_shell_sample"] = dist
        tmp["is_self"] = tmp["sample_id"].astype(str).eq(sid)
        tmp = tmp.sort_values(["distance_to_shell_sample", "sample_id"]).reset_index(drop=True)

        neigh = tmp.loc[~tmp["is_self"]].head(args.k).copy()
        neigh["shell_sample_id"] = sid
        neigh["k"] = args.k
        neighbor_rows.append(
            neigh[
                [
                    "shell_sample_id",
                    "sample_id",
                    "distance_to_shell_sample",
                    "state_domain",
                    "placeholder_condition_cur",
                    "biological_condition",
                    "tumor_status",
                    axis_map["H"],
                    axis_map["S"],
                    axis_map["M"],
                    axis_map["R"],
                    axis_map["phi"],
                ]
            ].rename(
                columns={
                    axis_map["H"]: "H_source",
                    axis_map["S"]: "S_source",
                    axis_map["M"]: "M_source",
                    axis_map["R"]: "R_source",
                    axis_map["phi"]: "phi_source",
                }
            )
        )

        counts_domain = neigh["state_domain"].astype(str).value_counts(dropna=False).to_dict()
        counts_cond = neigh["placeholder_condition_cur"].astype(str).value_counts(dropna=False).to_dict()

        summary_rows.append(
            {
                "shell_sample_id": sid,
                "k": args.k,
                "neighbor_median_distance": median_or_nan(neigh["distance_to_shell_sample"]),
                "neighbor_phi_median": median_or_nan(neigh[axis_map["phi"]]),
                "neighbor_H_median": median_or_nan(neigh[axis_map["H"]]),
                "neighbor_S_median": median_or_nan(neigh[axis_map["S"]]),
                "neighbor_M_median": median_or_nan(neigh[axis_map["M"]]),
                "neighbor_R_median": median_or_nan(neigh[axis_map["R"]]),
                "neighbor_state_domain_counts_json": json.dumps(counts_domain, sort_keys=True),
                "neighbor_condition_counts_json": json.dumps(counts_cond, sort_keys=True),
            }
        )

    neighbors_out = pd.concat(neighbor_rows, ignore_index=True)
    neighbors_out.to_csv(out_dir / "stable_shell_k_nearest_neighbors.csv", index=False)

    summary_out = pd.DataFrame(summary_rows).sort_values("shell_sample_id").reset_index(drop=True)
    summary_out.to_csv(out_dir / "stable_shell_neighbor_summary.csv", index=False)

    pooled = neighbors_out.copy()
    pooled["is_mgmt_methylated_neighbor"] = pooled["placeholder_condition_cur"].astype(str).eq("mgmt_methylated")
    pooled_counts = (
        pooled.groupby(["shell_sample_id", "state_domain"], dropna=False)
        .size()
        .reset_index(name="n_neighbors")
        .sort_values(["shell_sample_id", "n_neighbors"], ascending=[True, False])
        .reset_index(drop=True)
    )
    pooled_counts.to_csv(out_dir / "stable_shell_neighbor_domain_counts.csv", index=False)

    mgmt_counts = (
        pooled.groupby(["shell_sample_id", "placeholder_condition_cur"], dropna=False)
        .size()
        .reset_index(name="n_neighbors")
        .sort_values(["shell_sample_id", "n_neighbors"], ascending=[True, False])
        .reset_index(drop=True)
    )
    mgmt_counts.to_csv(out_dir / "stable_shell_neighbor_condition_counts.csv", index=False)

    meta = {
        "k": args.k,
        "n_total_samples_with_complete_axes": int(len(work)),
        "n_shell_samples_found": int(len(shell_work)),
        "shell_sample_ids": sorted(shell_work["sample_id"].astype(str).tolist()),
        "axis_column_map": axis_map,
    }
    with open(out_dir / "stable_shell_neighborhood_summary.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    print(f"[ok] wrote stable-shell neighborhood outputs to {out_dir}")
    print("[info] summary:")
    print(f"  k: {args.k}")
    print(f"  n_shell_samples_found: {len(shell_work)}")
    print(f"  shell_sample_ids: {sorted(shell_work['sample_id'].astype(str).tolist())}")


if __name__ == "__main__":
    main()
