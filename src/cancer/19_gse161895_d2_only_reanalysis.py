#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml
from scipy.stats import chi2_contingency, mannwhitneyu
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


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


def mwu_test(a: pd.Series, b: pd.Series) -> dict:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) == 0 or len(b) == 0:
        return {
            "n_a": len(a),
            "n_b": len(b),
            "median_a": None,
            "median_b": None,
            "delta_median_a_minus_b": None,
            "u_statistic": None,
            "p_value": None,
        }
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(a.median()),
        "median_b": float(b.median()),
        "delta_median_a_minus_b": float(a.median() - b.median()),
        "u_statistic": float(u),
        "p_value": float(p),
    }


def main() -> None:
    cfg = read_config()
    proc_dir = REPO_ROOT / cfg["processed_dir"]
    state_dir = REPO_ROOT / cfg["results_dir"] / "state_domains"
    out_dir = REPO_ROOT / cfg["results_dir"] / "d2_only_reanalysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    state = pd.read_parquet(proc_dir / "state_table.parquet")
    domains = pd.read_parquet(state_dir / "state_domains.parquet")[["cell_id", "state_domain"]]
    meta = pd.read_parquet(proc_dir / "cell_metadata_registry_aligned.parquet")

    df = state.merge(domains, on="cell_id", how="inner").merge(
        meta[["cell_id", "patient_or_donor_id", "source_class", "treatment_state", "source_name_ch1", "characteristics_ch1"]],
        on="cell_id",
        how="left",
    )

    d2 = df[df["state_domain"] == "D2"].copy()

    X = d2[["H", "S", "M", "R"]].values

    sil_rows = []
    best_k = None
    best_score = -1
    best_labels = None

    for k in range(2, 7):
        km = KMeans(n_clusters=k, random_state=cfg["seed"], n_init=20)
        labels = km.fit_predict(X)
        score = silhouette_score(X, labels)
        sil_rows.append({"k": k, "silhouette": float(score)})
        if score > best_score:
            best_score = score
            best_k = k
            best_labels = labels

    sil_df = pd.DataFrame(sil_rows)
    sil_df.to_csv(out_dir / "d2_silhouette_scan.csv", index=False)

    d2["d2_subdomain"] = [f"D2_{int(x)+1}" for x in best_labels]
    d2.to_parquet(out_dir / "d2_state_domains.parquet", index=False)
    d2["d2_subdomain"].value_counts().rename_axis("d2_subdomain").reset_index(name="n").to_csv(
        out_dir / "d2_state_domain_counts.csv", index=False
    )

    # metadata enrichment inside D2
    enrichment_rows = []
    candidate_cols = ["patient_or_donor_id", "source_class", "treatment_state", "source_name_ch1", "characteristics_ch1"]
    for col in candidate_cols:
        tmp = d2[[col, "d2_subdomain"]].copy()
        tmp[col] = tmp[col].fillna("NA").astype(str)
        vc = tmp[col].value_counts()
        keep = vc.index[:20]
        tmp.loc[~tmp[col].isin(keep), col] = "OTHER"
        table = pd.crosstab(tmp[col], tmp["d2_subdomain"])
        table.to_csv(out_dir / f"{col}_by_d2_subdomain.csv")
        stats = safe_chi2(table)
        enrichment_rows.append({
            "metadata_column": col,
            "n_categories": int(table.shape[0]),
            "chi2": stats["chi2"],
            "p_value": stats["p_value"],
            "dof": stats["dof"],
        })

    enrich_df = pd.DataFrame(enrichment_rows).sort_values("p_value", na_position="last")
    enrich_df.to_csv(out_dir / "d2_metadata_enrichment_summary.csv", index=False)

    # if best_k == 2, do direct HRSM tests between the two D2 subdomains
    hrsm_rows = []
    unique_sub = sorted(d2["d2_subdomain"].dropna().unique().tolist())
    if len(unique_sub) == 2:
        a_label, b_label = unique_sub
        for metric in ["H", "S", "M", "R", "phi"]:
            a = d2.loc[d2["d2_subdomain"] == a_label, metric]
            b = d2.loc[d2["d2_subdomain"] == b_label, metric]
            res = mwu_test(a, b)
            res["metric"] = metric
            res["group_a"] = a_label
            res["group_b"] = b_label
            hrsm_rows.append(res)

    hrsm_df = pd.DataFrame(hrsm_rows)
    if not hrsm_df.empty:
        hrsm_df.to_csv(out_dir / "d2_subdomain_hrsm_tests.csv", index=False)

    summary = {
        "n_d2_cells": int(len(d2)),
        "best_k": int(best_k),
        "best_silhouette": float(best_score),
        "d2_subdomains": d2["d2_subdomain"].value_counts().to_dict(),
    }
    with open(out_dir / "d2_reanalysis_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))
    print("\n[D2 metadata enrichment]")
    print(enrich_df)
    if not hrsm_df.empty:
        print("\n[D2 subdomain HRSM tests]")
        print(hrsm_df)


if __name__ == "__main__":
    main()
