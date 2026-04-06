#!/usr/bin/env python3

from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

REPO = Path.home() / "cancer_cells_hrsm_sims"
DATASET = "gse161895"


def load_inputs():
    state = pd.read_parquet(REPO / "data" / "processed" / DATASET / "state_table.parquet")
    domains = pd.read_parquet(REPO / "results" / DATASET / "state_domains" / "state_domains.parquet")
    meta = pd.read_parquet(REPO / "data" / "metadata" / DATASET / "cell_metadata_registry.parquet")
    return state, domains, meta


def attach_metadata(state: pd.DataFrame, domains: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    state = state.copy()
    domains = domains.copy()
    meta = meta.copy()

    if "cell_id" not in state.columns:
        raise KeyError("state_table.parquet is missing cell_id")
    if "cell_id" not in domains.columns:
        raise KeyError("state_domains.parquet is missing cell_id")
    if "state_domain" not in domains.columns:
        raise KeyError("state_domains.parquet is missing state_domain")
    if "gsm" not in meta.columns:
        raise KeyError("cell_metadata_registry.parquet is missing gsm")

    state["cell_id"] = state["cell_id"].astype(str)
    domains["cell_id"] = domains["cell_id"].astype(str)
    meta["gsm"] = meta["gsm"].astype(str)

    state["gsm"] = state["cell_id"].str.split("_", n=1).str[0]

    merged = (
        state.merge(domains[["cell_id", "state_domain"]], on="cell_id", how="left")
             .merge(meta, on="gsm", how="left", suffixes=("", "_meta"))
    )

    debug = {
        "n_state_rows": int(len(state)),
        "n_domain_rows": int(len(domains)),
        "n_meta_rows": int(len(meta)),
        "n_rows_with_metadata": int(
            merged[["patient_or_donor_id", "source_class", "treatment_state"]]
            .notna().any(axis=1).sum()
        ),
        "n_rows_without_metadata": int(
            (~merged[["patient_or_donor_id", "source_class", "treatment_state"]]
             .notna().any(axis=1)).sum()
        ),
        "source_class_unique_head": sorted(
            merged["source_class"].fillna("").astype(str).unique().tolist()
        )[:20],
        "patient_or_donor_unique_head": sorted(
            merged["patient_or_donor_id"].fillna("").astype(str).unique().tolist()
        )[:20],
    }
    print(json.dumps(debug, indent=2))
    return merged


def donor_mask(df: pd.DataFrame) -> pd.Series:
    source_class = df["source_class"].fillna("").astype(str).str.lower()
    patient_label = df["patient_or_donor_id"].fillna("").astype(str).str.lower()
    source_name = df["source_name_ch1"].fillna("").astype(str).str.lower()

    mask = source_class.eq("normal_donor")
    mask = mask | patient_label.str.startswith("normaldonor")
    mask = mask | source_name.str.contains("normal", regex=False)
    return mask


def add_reversion_diagnostics(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    donor = df[donor_mask(df)].copy()
    if donor.empty:
        raise ValueError("No donor-like cells found after GSM-based metadata merge.")

    hrsm_cols = ["H", "S", "M", "R", "phi"]
    donor_centroid = donor[hrsm_cols].median()

    # Euclidean distance to donor centroid in HRSM space
    delta = df[hrsm_cols].subtract(donor_centroid, axis=1)
    df["reversion_distance_to_donor"] = np.sqrt((delta ** 2).sum(axis=1))

    # Plasticity proxy, distance away from own domain centroid
    domain_centroids = df.groupby("state_domain")[hrsm_cols].median()
    plasticity_vals = []
    for _, row in df.iterrows():
        centroid = domain_centroids.loc[row["state_domain"]]
        d = row[hrsm_cols] - centroid
        plasticity_vals.append(float(np.sqrt((d ** 2).sum())))
    df["plasticity_distance_to_domain"] = plasticity_vals

    return df


def summarize_tests(df: pd.DataFrame, metric: str, group_col: str, a: str, b: str) -> dict:
    xa = df.loc[df[group_col] == a, metric].dropna().to_numpy()
    xb = df.loc[df[group_col] == b, metric].dropna().to_numpy()
    u, p = mannwhitneyu(xa, xb, alternative="two-sided")
    return {
        "metric": metric,
        "group_a": a,
        "group_b": b,
        "n_a": int(len(xa)),
        "n_b": int(len(xb)),
        "median_a": float(np.median(xa)),
        "median_b": float(np.median(xb)),
        "delta_median_a_minus_b": float(np.median(xa) - np.median(xb)),
        "u_statistic": float(u),
        "p_value": float(p),
    }


def main():
    outdir = REPO / "results" / DATASET / "plasticity_reversion_diagnostics"
    outdir.mkdir(parents=True, exist_ok=True)

    state, domains, meta = load_inputs()
    df = attach_metadata(state, domains, meta)
    df = add_reversion_diagnostics(df)

    df.to_csv(outdir / "plasticity_reversion_per_cell.csv", index=False)

    tests = []
    for metric in ["plasticity_distance_to_domain", "reversion_distance_to_donor"]:
        tests.append(summarize_tests(df, metric, "state_domain", "D1", "D2"))

    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(outdir / "plasticity_reversion_domain_tests.csv", index=False)

    summary = {
        "n_cells": int(len(df)),
        "n_donor_like_cells": int(donor_mask(df).sum()),
        "outputs": {
            "per_cell_csv": str(outdir / "plasticity_reversion_per_cell.csv"),
            "tests_csv": str(outdir / "plasticity_reversion_domain_tests.csv"),
        },
    }
    with open(outdir / "plasticity_reversion_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(tests_df.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
