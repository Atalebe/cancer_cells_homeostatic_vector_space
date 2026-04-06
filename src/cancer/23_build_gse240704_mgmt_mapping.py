from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


PATTERN = re.compile(
    r"mgmt promoter methylation status\s*\[0=non-methylated,\s*1=methylated\]\s*:\s*([01])",
    flags=re.IGNORECASE,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Parse MGMT promoter methylation status from GSE240704 series-matrix metadata."
    )
    parser.add_argument("--series-metadata-csv", required=True)
    parser.add_argument("--mapping-template-csv", required=True)
    parser.add_argument("--outcsv", required=True)
    return parser.parse_args()


def parse_mgmt(value: object) -> tuple[str | None, str | None]:
    if pd.isna(value):
        return None, None
    text = str(value).strip()
    m = PATTERN.search(text)
    if not m:
        return None, None
    raw = m.group(1)
    if raw == "0":
        return "mgmt_non_methylated", "non_methylated"
    if raw == "1":
        return "mgmt_methylated", "methylated"
    return None, None


def main() -> None:
    args = parse_args()

    sm = pd.read_csv(args.series_metadata_csv, dtype=str).fillna("")
    mt = pd.read_csv(args.mapping_template_csv, dtype=str).fillna("")

    if "sample_id" not in sm.columns or "sample_id" not in mt.columns:
        raise ValueError("Both files must contain sample_id")

    mgmt_df = sm[["sample_id", "characteristics_ch1", "title", "geo_accession", "source_name_ch1"]].copy()

    parsed = mgmt_df["characteristics_ch1"].apply(parse_mgmt)
    mgmt_df["placeholder_condition"] = parsed.apply(lambda x: x[0] if x else None)
    mgmt_df["tumor_status"] = parsed.apply(lambda x: x[1] if x else None)

    mgmt_df["placeholder_group"] = "glioblastoma"
    mgmt_df["disease_group"] = "glioblastoma"
    mgmt_df["cell_type"] = mgmt_df["source_name_ch1"].replace("", "Glioblastoma")
    mgmt_df["placeholder_note"] = mgmt_df["characteristics_ch1"]

    merged = mt.drop(
        columns=[c for c in [
            "placeholder_condition",
            "placeholder_group",
            "placeholder_note",
            "disease_group",
            "cell_type",
            "tumor_status",
        ] if c in mt.columns],
        errors="ignore",
    ).merge(
        mgmt_df[[
            "sample_id",
            "placeholder_condition",
            "placeholder_group",
            "placeholder_note",
            "disease_group",
            "cell_type",
            "tumor_status",
            "characteristics_ch1",
            "title",
            "geo_accession",
        ]],
        on="sample_id",
        how="left",
    )

    # Keep a simple biological label too
    merged["biological_condition"] = merged["placeholder_condition"]

    # Natural sort by sample number if present
    if "sample_number" in merged.columns:
        merged["__sample_number__"] = pd.to_numeric(merged["sample_number"], errors="coerce")
        merged = merged.sort_values(["__sample_number__", "sample_id"]).drop(columns="__sample_number__")
    else:
        merged = merged.sort_values("sample_id")

    Path(args.outcsv).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.outcsv, index=False)

    vc = merged["placeholder_condition"].fillna("missing").value_counts(dropna=False).to_dict()

    print("[ok] wrote MGMT-curated mapping:", args.outcsv)
    print("[info] rows:", len(merged))
    print("[info] condition counts:")
    for k, v in vc.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
