#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml
from scipy.stats import chi2_contingency


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def safe_chi2(table: pd.DataFrame) -> dict:
    if table.shape[0] < 2 or table.shape[1] < 2:
        return {"chi2": None, "p_value": None, "dof": None}
    try:
        chi2, p, dof, _ = chi2_contingency(table)
        return {"chi2": float(chi2), "p_value": float(p), "dof": int(dof)}
    except Exception:
        return {"chi2": None, "p_value": None, "dof": None}


def main() -> None:
    cfg = read_config()
    proc_dir = REPO_ROOT / cfg["processed_dir"]
    state_dir = REPO_ROOT / cfg["results_dir"] / "state_domains"
    out_dir = REPO_ROOT / cfg["results_dir"] / "domain_metadata_enrichment"
    out_dir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_parquet(proc_dir / "cell_metadata_registry_aligned.parquet")
    states = pd.read_parquet(state_dir / "state_domains.parquet")[["cell_id", "state_domain"]]

    df = meta.merge(states, on="cell_id", how="inner")

    candidate_cols = [
        "source_class",
        "treatment_state",
        "patient_or_donor_id",
        "source_name_ch1",
        "characteristics_ch1",
        "title",
    ]

    test_rows = []

    for col in candidate_cols:
        if col not in df.columns:
            continue

        tmp = df[[col, "state_domain"]].copy()
        tmp[col] = tmp[col].fillna("NA").astype(str)

        # keep top categories only for huge text-heavy cols
        vc = tmp[col].value_counts()
        keep = vc.index[:20]
        tmp.loc[~tmp[col].isin(keep), col] = "OTHER"

        table = pd.crosstab(tmp[col], tmp["state_domain"])
        table.to_csv(out_dir / f"{col}_by_state_domain.csv")

        stats = safe_chi2(table)
        test_rows.append(
            {
                "metadata_column": col,
                "n_categories": int(table.shape[0]),
                "chi2": stats["chi2"],
                "p_value": stats["p_value"],
                "dof": stats["dof"],
            }
        )

    summary_df = pd.DataFrame(test_rows).sort_values("p_value", na_position="last")
    summary_df.to_csv(out_dir / "domain_metadata_enrichment_summary.csv", index=False)

    with open(out_dir / "domain_metadata_enrichment_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary_df.to_dict(orient="records"), f, indent=2)

    print(summary_df)


if __name__ == "__main__":
    main()
