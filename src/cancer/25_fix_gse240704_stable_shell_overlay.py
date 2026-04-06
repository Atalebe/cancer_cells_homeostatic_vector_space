from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DOMAIN_CANDIDATES = [
    "results/gse240704/state_domains/state_domain_assignments.csv",
    "results/gse240704/state_domains/sample_state_domains.csv",
    "results/gse240704/state_domains/state_domains.csv",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Repair stable-shell overlay by attaching state-domain labels and auditing missing annotations."
    )
    parser.add_argument(
        "--overlay-csv",
        required=True,
        help="Path to stable_shell_mgmt_overlay.csv",
    )
    parser.add_argument(
        "--stable-shell-csv",
        required=True,
        help="Path to stable_shell_samples_inspection.csv",
    )
    parser.add_argument(
        "--sample-annotations-curated",
        required=True,
        help="Curated sample annotation parquet/csv",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory",
    )
    return parser.parse_args()


def load_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(path)
    if p.suffix.lower() == ".parquet":
        return pd.read_parquet(p)
    return pd.read_csv(p)


def load_domain_table() -> tuple[pd.DataFrame, str]:
    for p in DOMAIN_CANDIDATES:
        path = Path(p)
        if not path.exists():
            continue
        df = pd.read_csv(path)
        cols = set(df.columns)
        if {"sample_id", "state_domain"}.issubset(cols):
            return df[["sample_id", "state_domain"]].drop_duplicates("sample_id"), p
    raise FileNotFoundError(
        "No usable state-domain assignment file found in expected locations."
    )


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    overlay = pd.read_csv(args.overlay_csv)
    stable = pd.read_csv(args.stable_shell_csv)
    ann = load_table(args.sample_annotations_curated)
    dom, dom_path = load_domain_table()

    repaired = (
        stable[["sample_id", "phi", "H", "S", "M", "R"]]
        .drop_duplicates("sample_id")
        .merge(dom, on="sample_id", how="left")
        .merge(
            ann.drop_duplicates("sample_id"),
            on="sample_id",
            how="left",
            suffixes=("", "_ann"),
        )
    )

    keep_cols = [
        "sample_id",
        "state_domain",
        "placeholder_condition_cur",
        "biological_condition",
        "tumor_status",
        "characteristics_ch1",
        "phi",
        "H",
        "S",
        "M",
        "R",
    ]
    keep_cols = [c for c in keep_cols if c in repaired.columns]
    repaired = repaired[keep_cols]

    repaired.to_csv(outdir / "stable_shell_mgmt_overlay_repaired.csv", index=False)

    missing = repaired[repaired.isna().any(axis=1)].copy()
    missing.to_csv(outdir / "stable_shell_missing_annotation_audit.csv", index=False)

    pd.DataFrame(
        {
            "domain_source_file": [dom_path],
            "n_overlay_rows_input": [len(overlay)],
            "n_overlay_rows_repaired": [len(repaired)],
            "n_rows_with_any_missing_fields": [len(missing)],
            "n_missing_condition": [
                int(repaired["placeholder_condition_cur"].isna().sum())
                if "placeholder_condition_cur" in repaired.columns
                else -1
            ],
            "n_missing_state_domain": [
                int(repaired["state_domain"].isna().sum())
                if "state_domain" in repaired.columns
                else -1
            ],
        }
    ).to_csv(outdir / "stable_shell_overlay_repair_summary.csv", index=False)

    print("[ok] wrote repaired overlay:", outdir / "stable_shell_mgmt_overlay_repaired.csv")
    print("[ok] wrote missing annotation audit:", outdir / "stable_shell_missing_annotation_audit.csv")
    print("[ok] wrote repair summary:", outdir / "stable_shell_overlay_repair_summary.csv")
    print("[info] state-domain source:", dom_path)


if __name__ == "__main__":
    main()
