from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inspect stable-shell and irreversible-shell samples for GSE240704."
    )
    parser.add_argument(
        "--state-table",
        required=True,
        help="Path to sample-level HRSM state table parquet.",
    )
    parser.add_argument(
        "--sample-annotations",
        required=True,
        help="Path to merged sample annotations parquet.",
    )
    parser.add_argument(
        "--candidate-shell-table",
        required=True,
        help="Path to candidate shell table parquet or csv.",
    )
    parser.add_argument(
        "--state-domains",
        required=True,
        help="Path to state domains parquet or csv.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for inspection tables.",
    )
    return parser.parse_args()


def read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() == ".parquet":
        return pd.read_parquet(p)
    if p.suffix.lower() == ".csv":
        return pd.read_csv(p)
    raise ValueError(f"Unsupported file format: {p}")


def ensure_sample_id(df: pd.DataFrame, name: str) -> pd.DataFrame:
    if "sample_id" not in df.columns:
        raise ValueError(f"{name} is missing required column: sample_id")
    return df.copy()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    state = ensure_sample_id(read_table(args.state_table), "state_table")
    ann = ensure_sample_id(read_table(args.sample_annotations), "sample_annotations")
    shell = ensure_sample_id(read_table(args.candidate_shell_table), "candidate_shell_table")
    dom = ensure_sample_id(read_table(args.state_domains), "state_domains")

    merged = (
        state.merge(ann, on="sample_id", how="left", suffixes=("", "_ann"))
             .merge(dom, on="sample_id", how="left", suffixes=("", "_dom"))
             .merge(shell, on="sample_id", how="left", suffixes=("", "_shell"))
    )

    stable_flag_candidates = [
        "is_candidate_stable_shell",
        "candidate_stable_shell",
        "stable_shell_flag",
    ]
    irreversible_flag_candidates = [
        "is_candidate_irreversible_shell",
        "candidate_irreversible_shell",
        "irreversible_shell_flag",
    ]

    stable_flag = next((c for c in stable_flag_candidates if c in merged.columns), None)
    irreversible_flag = next((c for c in irreversible_flag_candidates if c in merged.columns), None)

    if stable_flag is None:
        raise ValueError(
            "Could not find a stable shell flag column. "
            f"Checked: {stable_flag_candidates}"
        )

    stable_samples = merged.loc[merged[stable_flag].fillna(False).astype(bool)].copy()

    if irreversible_flag is not None:
        irreversible_samples = merged.loc[
            merged[irreversible_flag].fillna(False).astype(bool)
        ].copy()
    else:
        irreversible_samples = merged.iloc[0:0].copy()

    preferred_cols = [
        "sample_id",
        "phi",
        "H",
        "S",
        "M",
        "R",
        "state_domain",
        stable_flag,
    ]
    if irreversible_flag is not None:
        preferred_cols.append(irreversible_flag)

    for extra in [
        "sample_number",
        "placeholder_condition",
        "placeholder_group",
        "placeholder_batch",
        "placeholder_note",
    ]:
        if extra in stable_samples.columns:
            preferred_cols.append(extra)

    preferred_cols = [c for c in preferred_cols if c in stable_samples.columns]

    stable_samples = stable_samples.sort_values(["phi", "S", "M"], ascending=[False, False, False])
    stable_samples[preferred_cols].to_csv(
        outdir / "stable_shell_samples_inspection.csv",
        index=False,
    )

    irreversible_samples = irreversible_samples.sort_values(
        ["phi", "S", "M"], ascending=[False, False, False]
    )
    irr_cols = [c for c in preferred_cols if c in irreversible_samples.columns]
    irreversible_samples[irr_cols].to_csv(
        outdir / "irreversible_shell_samples_inspection.csv",
        index=False,
    )

    summary = {
        "n_total_samples": int(len(merged)),
        "n_stable_shell_samples": int(len(stable_samples)),
        "n_irreversible_shell_samples": int(len(irreversible_samples)),
        "stable_shell_flag_column": stable_flag,
        "irreversible_shell_flag_column": irreversible_flag,
        "stable_shell_sample_ids": stable_samples["sample_id"].tolist(),
        "irreversible_shell_sample_ids": irreversible_samples["sample_id"].tolist(),
    }

    with open(outdir / "stable_shell_inspection_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote stable-shell inspection tables to", outdir)
    print("[info] n_stable_shell_samples:", len(stable_samples))
    print("[info] n_irreversible_shell_samples:", len(irreversible_samples))
    if len(stable_samples) > 0:
        print("[info] stable shell sample IDs:")
        for sid in stable_samples["sample_id"].tolist():
            print(" ", sid)


if __name__ == "__main__":
    main()
