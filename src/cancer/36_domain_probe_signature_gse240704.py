#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Probe-level domain signature tests for D2 and D3.")
    p.add_argument("--selected-probe-matrix", required=True, help="Parquet with ID_REF + SAMPLE columns")
    p.add_argument("--merged-curated-csv", required=True, help="Merged curated domain-condition table")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--top-n", type=int, default=200, help="Top probes to export per contrast")
    return p.parse_args()


def safe_mwu(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) == 0 or len(b) == 0:
        return np.nan, np.nan
    try:
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        return float(u), float(p)
    except Exception:
        return np.nan, np.nan


def rank_biserial_from_u(u: float, n_a: int, n_b: int) -> float:
    if not np.isfinite(u) or n_a == 0 or n_b == 0:
        return np.nan
    return float((2.0 * u / (n_a * n_b)) - 1.0)


def benjamini_hochberg(pvals: pd.Series) -> pd.Series:
    s = pvals.astype(float).copy()
    out = pd.Series(np.nan, index=s.index, dtype=float)
    mask = s.notna()
    if mask.sum() == 0:
        return out
    pv = s[mask].values
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out.loc[mask] = q[np.argsort(order)]
    return out


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    merged = pd.read_csv(args.merged_curated_csv)
    mat = pd.read_parquet(args.selected_probe_matrix)

    if "ID_REF" not in mat.columns:
        raise SystemExit("selected probe matrix must contain ID_REF")

    sample_cols = [c for c in mat.columns if c != "ID_REF"]
    sample_set = set(sample_cols)

    if "sample_id" not in merged.columns:
        raise SystemExit("merged curated table must contain sample_id")
    if "state_domain" not in merged.columns:
        raise SystemExit("merged curated table must contain state_domain")

    merged = merged.copy()
    merged["sample_id"] = merged["sample_id"].astype(str)

    domain_samples = {
        d: [s for s in merged.loc[merged["state_domain"] == d, "sample_id"].tolist() if s in sample_set]
        for d in sorted(merged["state_domain"].dropna().unique())
    }

    contrasts = [
        ("D3", "not_D3", domain_samples.get("D3", []), [s for s in sample_cols if s not in set(domain_samples.get("D3", []))]),
        ("D2", "not_D2", domain_samples.get("D2", []), [s for s in sample_cols if s not in set(domain_samples.get("D2", []))]),
        ("D3", "D2", domain_samples.get("D3", []), domain_samples.get("D2", [])),
        ("D2", "D1", domain_samples.get("D2", []), domain_samples.get("D1", [])),
        ("D3", "D1", domain_samples.get("D3", []), domain_samples.get("D1", [])),
    ]

    all_rows = []
    top_rows = []

    for label_a, label_b, cols_a, cols_b in contrasts:
        if len(cols_a) == 0 or len(cols_b) == 0:
            continue

        block = mat[["ID_REF"] + cols_a + cols_b].copy()
        arr_a = block[cols_a].to_numpy(dtype=float)
        arr_b = block[cols_b].to_numpy(dtype=float)

        med_a = np.nanmedian(arr_a, axis=1)
        med_b = np.nanmedian(arr_b, axis=1)
        mean_a = np.nanmean(arr_a, axis=1)
        mean_b = np.nanmean(arr_b, axis=1)
        std_a = np.nanstd(arr_a, axis=1)
        std_b = np.nanstd(arr_b, axis=1)

        rows = []
        for i, probe in enumerate(block["ID_REF"].astype(str).tolist()):
            u, p = safe_mwu(arr_a[i, :], arr_b[i, :])
            rb = rank_biserial_from_u(u, len(cols_a), len(cols_b))
            rows.append(
                {
                    "contrast_a": label_a,
                    "contrast_b": label_b,
                    "ID_REF": probe,
                    "n_a": len(cols_a),
                    "n_b": len(cols_b),
                    "median_a": float(med_a[i]),
                    "median_b": float(med_b[i]),
                    "delta_median_a_minus_b": float(med_a[i] - med_b[i]),
                    "mean_a": float(mean_a[i]),
                    "mean_b": float(mean_b[i]),
                    "delta_mean_a_minus_b": float(mean_a[i] - mean_b[i]),
                    "std_a": float(std_a[i]),
                    "std_b": float(std_b[i]),
                    "u_statistic": u,
                    "p_value": p,
                    "rank_biserial": rb,
                }
            )

        df = pd.DataFrame(rows)
        df["q_value_bh"] = benjamini_hochberg(df["p_value"])
        df["abs_delta_median"] = df["delta_median_a_minus_b"].abs()
        df["abs_rank_biserial"] = df["rank_biserial"].abs()

        all_rows.append(df)

        top_df = (
            df.sort_values(["q_value_bh", "abs_delta_median", "abs_rank_biserial"], ascending=[True, False, False])
              .head(args.top_n)
              .reset_index(drop=True)
        )
        top_rows.append(top_df)

        contrast_slug = f"{label_a}_vs_{label_b}"
        df.to_csv(outdir / f"{contrast_slug}_all_probes.csv", index=False)
        top_df.to_csv(outdir / f"{contrast_slug}_top_{args.top_n}_probes.csv", index=False)

    if not all_rows:
        raise SystemExit("no valid contrasts could be computed")

    all_df = pd.concat(all_rows, ignore_index=True)
    top_df = pd.concat(top_rows, ignore_index=True)

    all_df.to_csv(outdir / "domain_probe_signature_tests_all.csv", index=False)
    top_df.to_csv(outdir / "domain_probe_signature_tests_top.csv", index=False)

    summary = {
        "domains_present": sorted(merged["state_domain"].dropna().unique().tolist()),
        "contrast_count": int(len(top_rows)),
        "selected_probe_count": int(len(mat)),
        "sample_count": int(len(sample_cols)),
        "domain_sample_counts": {k: int(len(v)) for k, v in domain_samples.items()},
        "focus": "probe-level separation for D2 and D3, including D3 vs rest and D3 vs D2",
    }
    with open(outdir / "domain_probe_signature_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("[ok] wrote domain probe signature outputs to", outdir)
    print("[info] summary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
