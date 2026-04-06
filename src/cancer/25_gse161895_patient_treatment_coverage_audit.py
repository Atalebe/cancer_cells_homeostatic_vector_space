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


def chi2_from_crosstab(df: pd.DataFrame) -> dict:
    if df.shape[0] < 2 or df.shape[1] < 2:
        return {"chi2": None, "p_value": None, "dof": None}
    chi2, p, dof, _ = chi2_contingency(df)
    return {"chi2": float(chi2), "p_value": float(p), "dof": int(dof)}


def resolve_state_domain_column(df: pd.DataFrame) -> str:
    candidates = ["state_domain", "domain", "cluster", "label"]
    for col in candidates:
        if col in df.columns:
            return col
    raise ValueError(f"Could not find state-domain column in columns: {list(df.columns)}")


def main() -> None:
    cfg = read_config()

    proc_dir = REPO_ROOT / cfg["processed_dir"]
    results_dir = REPO_ROOT / cfg["results_dir"]
    d2_dir = results_dir / "d2_only_reanalysis"
    state_domains_path = results_dir / "state_domains" / "state_domains.parquet"

    out_dir = results_dir / "patient_treatment_coverage_audit"
    out_dir.mkdir(parents=True, exist_ok=True)

    state = pd.read_parquet(proc_dir / "state_table.parquet")[["cell_id", "H", "S", "M", "R", "phi"]]

    state_domains = pd.read_parquet(state_domains_path)
    state_domain_col = resolve_state_domain_column(state_domains)
    state_domains = state_domains.rename(columns={state_domain_col: "state_domain"})
    state_domains = state_domains[["cell_id", "state_domain"]]

    meta = pd.read_parquet(proc_dir / "cell_metadata_registry_aligned.parquet")[
        ["cell_id", "patient_or_donor_id", "source_class", "treatment_state", "source_name_ch1", "characteristics_ch1"]
    ]

    d2 = pd.read_parquet(d2_dir / "d2_state_domains.parquet")[["cell_id", "d2_subdomain"]]

    merged = (
        state
        .merge(state_domains, on="cell_id", how="left")
        .merge(meta, on="cell_id", how="left")
        .merge(d2, on="cell_id", how="left")
    )

    merged["treatment_state_clean"] = merged["treatment_state"].fillna("missing")

    patient_treatment = pd.crosstab(merged["patient_or_donor_id"], merged["treatment_state_clean"])
    patient_treatment.to_csv(out_dir / "patient_by_treatment_state.csv")

    sourceclass_treatment = pd.crosstab(merged["source_class"], merged["treatment_state_clean"])
    sourceclass_treatment.to_csv(out_dir / "source_class_by_treatment_state.csv")

    top_state_treatment = pd.crosstab(merged["state_domain"], merged["treatment_state_clean"])
    top_state_treatment.to_csv(out_dir / "state_domain_by_treatment_state.csv")

    d2_only = merged.loc[merged["d2_subdomain"].notna()].copy()
    d2_treatment = pd.crosstab(d2_only["d2_subdomain"], d2_only["treatment_state_clean"])
    d2_treatment.to_csv(out_dir / "d2_subdomain_by_treatment_state.csv")

    d2_1 = merged.loc[merged["d2_subdomain"] == "D2_1"].copy()
    d2_1_patient_treatment = pd.crosstab(d2_1["patient_or_donor_id"], d2_1["treatment_state_clean"])
    d2_1_patient_treatment.to_csv(out_dir / "d2_1_patient_by_treatment_state.csv")

    no_treated_entities = []
    if "treated" in patient_treatment.columns:
        for entity, row in patient_treatment.iterrows():
            if row.get("treated", 0) == 0:
                no_treated_entities.append(str(entity))
    else:
        no_treated_entities = patient_treatment.index.astype(str).tolist()

    only_untreated_entities = []
    for entity, row in patient_treatment.iterrows():
        treated = row.get("treated", 0)
        untreated = row.get("untreated", 0)
        missing = row.get("missing", 0)
        if treated == 0 and untreated > 0 and missing == 0:
            only_untreated_entities.append(str(entity))

    summary = {
        "n_cells_total": int(len(merged)),
        "n_unique_patients_or_donors": int(merged["patient_or_donor_id"].nunique(dropna=True)),
        "state_domain_source": str(state_domains_path),
        "no_treated_entities": no_treated_entities,
        "only_untreated_entities": only_untreated_entities,
        "patient_treatment_chi2": chi2_from_crosstab(patient_treatment),
        "sourceclass_treatment_chi2": chi2_from_crosstab(sourceclass_treatment),
        "state_domain_treatment_chi2": chi2_from_crosstab(top_state_treatment),
        "d2_subdomain_treatment_chi2": chi2_from_crosstab(d2_treatment),
        "outputs": {
            "patient_by_treatment": str(out_dir / "patient_by_treatment_state.csv"),
            "source_class_by_treatment": str(out_dir / "source_class_by_treatment_state.csv"),
            "state_domain_by_treatment": str(out_dir / "state_domain_by_treatment_state.csv"),
            "d2_subdomain_by_treatment": str(out_dir / "d2_subdomain_by_treatment_state.csv"),
            "d2_1_patient_by_treatment": str(out_dir / "d2_1_patient_by_treatment_state.csv"),
        },
    }

    with open(out_dir / "patient_treatment_coverage_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))
    print("\n[patient by treatment]")
    print(patient_treatment)
    print("\n[d2_1 patient by treatment]")
    print(d2_1_patient_treatment)


if __name__ == "__main__":
    main()
