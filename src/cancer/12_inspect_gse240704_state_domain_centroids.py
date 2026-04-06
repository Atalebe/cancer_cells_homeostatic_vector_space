from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


STATE_COLS = ["H", "S", "M", "R", "phi"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inspect centroid signatures of GSE240704 state domains."
    )
    parser.add_argument(
        "--state-table",
        required=True,
        help="Path to sample-level HRSM state table parquet.",
    )
    parser.add_argument(
        "--state-domains",
        required=True,
        help="Path to state domains parquet or csv.",
    )
    parser.add_argument(
        "--sample-annotations",
        required=False,
        default=None,
        help="Optional merged sample annotations parquet or csv.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for centroid inspection tables.",
    )
    return parser.parse_args()


def read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() == ".parquet":
        return pd.read_parquet(p)
    if p.suffix.lower() == ".csv":
        return pd.read_csv(p)
    raise ValueError(f"Unsupported file format: {p}")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    state = read_table(args.state_table)
    dom = read_table(args.state_domains)

    if "sample_id" not in state.columns or "sample_id" not in dom.columns:
        raise ValueError("Both state table and state domains must contain sample_id.")

    if "state_domain" not in dom.columns:
        raise ValueError("State domains table must contain state_domain.")

    df = state.merge(dom[["sample_id", "state_domain"]], on="sample_id", how="inner")

    missing = [c for c in STATE_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing state columns: {missing}")

    centroid = (
        df.groupby("state_domain")[STATE_COLS]
        .agg(["mean", "median", "std", "min", "max"])
        .round(6)
    )
    centroid.to_csv(outdir / "state_domain_centroid_signatures.csv")

    counts = (
        df.groupby("state_domain")
        .size()
        .rename("n")
        .reset_index()
        .sort_values("n", ascending=False)
    )
    counts.to_csv(outdir / "state_domain_counts_recomputed.csv", index=False)

    ranked = df.copy()
    ranked["abs_H"] = ranked["H"].abs()
    ranked["abs_S"] = ranked["S"].abs()
    ranked["abs_M"] = ranked["M"].abs()
    ranked["abs_R"] = ranked["R"].abs()

    domain_representatives = []
    for domain, sub in ranked.groupby("state_domain"):
        top_phi = sub.sort_values("phi", ascending=False).head(5).copy()
        top_phi.insert(0, "representative_type", "top_phi")
        low_phi = sub.sort_values("phi", ascending=True).head(5).copy()
        low_phi.insert(0, "representative_type", "low_phi")
        domain_representatives.append(top_phi)
        domain_representatives.append(low_phi)

    reps = pd.concat(domain_representatives, ignore_index=True)
    rep_cols = [c for c in ["representative_type", "state_domain", "sample_id", "phi", "H", "S", "M", "R"] if c in reps.columns]
    reps[rep_cols].to_csv(outdir / "state_domain_representative_samples.csv", index=False)

    signature_rows = []
    for domain, sub in df.groupby("state_domain"):
        row = {"state_domain": domain, "n": len(sub)}
        means = sub[STATE_COLS].mean()
        for col in STATE_COLS:
            row[f"{col}_mean"] = round(float(means[col]), 6)
        row["dominant_axis_by_abs_mean"] = max(
            ["H", "S", "M", "R"],
            key=lambda x: abs(row[f"{x}_mean"])
        )
        signature_rows.append(row)

    signature = pd.DataFrame(signature_rows).sort_values("n", ascending=False)
    signature.to_csv(outdir / "state_domain_signature_summary.csv", index=False)

    if args.sample_annotations:
        ann = read_table(args.sample_annotations)
        if "sample_id" in ann.columns:
            merged = df.merge(ann, on="sample_id", how="left")
            merged.to_parquet(outdir / "state_domain_state_plus_annotations.parquet", index=False)

    print("[ok] wrote centroid signature tables to", outdir)
    print("[info] domains:", ", ".join(sorted(df["state_domain"].dropna().astype(str).unique())))


if __name__ == "__main__":
    main()
