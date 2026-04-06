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

    state["cell_id"] = state["cell_id"].astype(str)
    domains["cell_id"] = domains["cell_id"].astype(str)
    meta["gsm"] = meta["gsm"].astype(str)

    state["gsm"] = state["cell_id"].str.split("_", n=1).str[0]

    df = (
        state.merge(domains[["cell_id", "state_domain"]], on="cell_id", how="left", validate="one_to_one")
             .merge(meta, on="gsm", how="left", validate="many_to_one")
    )

    df = df.loc[:, ~df.columns.duplicated()].copy()
    return df


def rank01(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    if s.nunique(dropna=True) <= 1:
        return pd.Series(np.full(len(s), 0.5), index=s.index)
    return s.rank(method="average", pct=True)


def add_state_age_proxy(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    h_rank = rank01(df["H"])
    m_rank = rank01(df["M"])
    phi_rank = rank01(df["phi"])
    r_rev_rank = 1.0 - rank01(df["R"])

    df["state_age_proxy"] = (h_rank + m_rank + phi_rank + r_rev_rank) / 4.0
    return df


def compare(df: pd.DataFrame, group_col: str, a: str, b: str, metric: str) -> dict:
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
    outdir = REPO / "results" / DATASET / "age_axis_proxy"
    outdir.mkdir(parents=True, exist_ok=True)

    state, domains, meta = load_inputs()
    df = attach_metadata(state, domains, meta)
    df = add_state_age_proxy(df)

    df.to_csv(outdir / "state_age_proxy_per_cell.csv", index=False)

    tests = [
        compare(df, "state_domain", "D1", "D2", "state_age_proxy"),
    ]

    d2_path = REPO / "results" / DATASET / "d2_only_reanalysis" / "d2_state_domains.parquet"
    if d2_path.exists():
        d2 = pd.read_parquet(d2_path).copy()

        if "cell_id" not in d2.columns:
            raise KeyError("d2_state_domains.parquet is missing cell_id")

        subdomain_col = None
        for candidate in ["d2_subdomain", "state_domain"]:
            if candidate in d2.columns:
                subdomain_col = candidate
                break
        if subdomain_col is None:
            raise KeyError("No d2 subdomain column found in d2_state_domains.parquet")

        d2 = d2[["cell_id", subdomain_col]].copy().rename(columns={subdomain_col: "d2_subdomain"})

        if "d2_subdomain" in df.columns:
            df = df.drop(columns=["d2_subdomain"])

        df2 = df.merge(d2, on="cell_id", how="left", validate="one_to_one")
        df2 = df2.loc[:, ~df2.columns.duplicated()].copy()

        sub = df2[df2["d2_subdomain"].astype(str) == "D2_1"].copy()
        sub["treatment_state"] = sub["treatment_state"].fillna("").astype(str)

        if {"treated", "untreated"}.issubset(set(sub["treatment_state"].unique())):
            tests.append(compare(sub, "treatment_state", "treated", "untreated", "state_age_proxy"))

    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(outdir / "state_age_proxy_tests.csv", index=False)

    summary = {
        "definition": "internal state-age proxy = mean(rank(H), rank(M), rank(phi), 1-rank(R))",
        "note": "This is not chronological patient age. It is a state-position age proxy.",
        "outputs": {
            "per_cell_csv": str(outdir / "state_age_proxy_per_cell.csv"),
            "tests_csv": str(outdir / "state_age_proxy_tests.csv"),
        },
    }
    with open(outdir / "state_age_proxy_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(tests_df.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
