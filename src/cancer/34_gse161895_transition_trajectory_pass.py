#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import mannwhitneyu


REPO = Path(__file__).resolve().parents[2]
RESULTS = REPO / "results" / "gse161895"
OUTDIR = RESULTS / "transition_trajectory_pass"
OUTDIR.mkdir(parents=True, exist_ok=True)


def load_state() -> pd.DataFrame:
    state = pd.read_parquet(REPO / "data" / "processed" / "gse161895" / "state_table.parquet")
    needed = ["cell_id", "H", "S", "M", "R", "phi"]
    missing = [c for c in needed if c not in state.columns]
    if missing:
        raise KeyError(f"Missing state columns: {missing}")
    return state[needed].copy()


def load_domains() -> pd.DataFrame:
    d1d2 = pd.read_parquet(RESULTS / "state_domains" / "state_domains.parquet")
    d2 = pd.read_parquet(RESULTS / "d2_only_reanalysis" / "d2_state_domains.parquet")
    keep_d1d2 = [c for c in ["cell_id", "state_domain"] if c in d1d2.columns]
    keep_d2 = [c for c in ["cell_id", "d2_subdomain"] if c in d2.columns]
    return d1d2[keep_d1d2].merge(d2[keep_d2], on="cell_id", how="left")


def load_meta() -> pd.DataFrame:
    meta = pd.read_parquet(REPO / "data" / "metadata" / "gse161895" / "cell_metadata_registry.parquet")
    state = load_state()[["cell_id"]].copy()
    gsm_map = state.copy()
    gsm_map["gsm"] = gsm_map["cell_id"].astype(str).str.split("_").str[0]
    key = "gsm" if "gsm" in meta.columns else ("geo_accession" if "geo_accession" in meta.columns else None)
    if key is None:
        raise KeyError("Metadata must contain gsm or geo_accession for GSE161895.")
    return gsm_map.merge(meta.rename(columns={key: "gsm"}), on="gsm", how="left").drop(columns=["gsm"])


def choose_group(row: pd.Series) -> str:
    source_class = str(row.get("source_class", ""))
    state_domain = str(row.get("state_domain", ""))
    d2_subdomain = str(row.get("d2_subdomain", ""))
    if source_class == "normal_donor":
        return "Donor"
    if state_domain == "D1":
        return "D1"
    if d2_subdomain == "D2_1":
        return "D2_1"
    if d2_subdomain == "D2_2":
        return "D2_2"
    if state_domain == "D2":
        return "D2_other"
    return "Other"


def project_to_axis(values: np.ndarray, start: np.ndarray, end: np.ndarray) -> np.ndarray:
    axis = end - start
    denom = np.dot(axis, axis)
    if denom == 0:
        return np.zeros(len(values))
    return np.dot(values - start, axis) / denom


def main() -> None:
    state = load_state()
    domains = load_domains()
    meta = load_meta()

    df = state.merge(domains, on="cell_id", how="left").merge(meta, on="cell_id", how="left")
    df["trajectory_group"] = df.apply(choose_group, axis=1)

    hrsm = ["H", "S", "M", "R"]
    usable_for_centroids = df[df["trajectory_group"].isin(["Donor", "D1", "D2_1", "D2_2"])].copy()
    centroids = usable_for_centroids.groupby("trajectory_group")[hrsm].median().loc[["Donor", "D2_2", "D2_1", "D1"]]

    dist_mat = cdist(df[hrsm].to_numpy(), centroids.to_numpy(), metric="euclidean")
    for i, grp in enumerate(centroids.index):
        df[f"dist_to_{grp}"] = dist_mat[:, i]
    df["nearest_trajectory_state"] = centroids.index.to_numpy()[dist_mat.argmin(axis=1)]

    donor = centroids.loc["Donor"].to_numpy()
    d1 = centroids.loc["D1"].to_numpy()
    d21 = centroids.loc["D2_1"].to_numpy()
    d22 = centroids.loc["D2_2"].to_numpy()
    vals = df[hrsm].to_numpy()

    df["trajectory_progress_donor_to_d1"] = project_to_axis(vals, donor, d1)
    df["trajectory_progress_donor_to_d2_1"] = project_to_axis(vals, donor, d21)
    df["trajectory_progress_d2_2_to_d2_1"] = project_to_axis(vals, d22, d21)
    df["d2_branch_bias"] = df["dist_to_D2_2"] - df["dist_to_D2_1"]
    df["donor_reversion_bias"] = df["dist_to_D1"] - df["dist_to_Donor"]

    usable = df[df["trajectory_group"].isin(["Donor", "D1", "D2_1", "D2_2"])].copy()
    tests = []
    for metric, ga, gb in [
        ("trajectory_progress_donor_to_d1", "D1", "D2_1"),
        ("trajectory_progress_donor_to_d2_1", "D2_1", "D2_2"),
        ("d2_branch_bias", "D2_1", "D2_2"),
        ("donor_reversion_bias", "D2_1", "D1"),
    ]:
        a = usable.loc[usable["trajectory_group"] == ga, metric].dropna()
        b = usable.loc[usable["trajectory_group"] == gb, metric].dropna()
        stat, p = mannwhitneyu(a, b, alternative="two-sided")
        tests.append({
            "metric": metric,
            "group_a": ga,
            "group_b": gb,
            "n_a": int(a.shape[0]),
            "n_b": int(b.shape[0]),
            "median_a": float(a.median()),
            "median_b": float(b.median()),
            "delta_median_a_minus_b": float(a.median() - b.median()),
            "u_statistic": float(stat),
            "p_value": float(p),
        })

    centroids.reset_index().to_csv(OUTDIR / "trajectory_centroids.csv", index=False)
    df.to_csv(OUTDIR / "trajectory_per_cell.csv", index=False)
    pd.DataFrame(tests).to_csv(OUTDIR / "trajectory_group_tests.csv", index=False)

    summary = {
        "trajectory_order": ["Donor", "D2_2", "D2_1", "D1"],
        "n_cells": int(df.shape[0]),
        "n_usable_cells": int(usable.shape[0]),
        "outputs": {
            "centroids_csv": str(OUTDIR / "trajectory_centroids.csv"),
            "per_cell_csv": str(OUTDIR / "trajectory_per_cell.csv"),
            "tests_csv": str(OUTDIR / "trajectory_group_tests.csv"),
        },
    }
    (OUTDIR / "transition_trajectory_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(pd.DataFrame(tests).to_string(index=False))


if __name__ == "__main__":
    main()
