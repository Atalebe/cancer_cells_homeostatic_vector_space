#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_jsonish_counts(s: str) -> dict:
    if pd.isna(s):
        return {}
    try:
        return json.loads(s)
    except Exception:
        return {}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--neighbor-summary-csv", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.neighbor_summary_csv)

    rows = []
    for _, r in df.iterrows():
        dom = parse_jsonish_counts(r.get("neighbor_state_domain_counts_json"))
        cond = parse_jsonish_counts(r.get("neighbor_condition_counts_json"))

        n_dom = sum(dom.values()) if len(dom) else 0
        n_cond = sum(cond.values()) if len(cond) else 0

        max_dom = max(dom.values()) if len(dom) else 0
        max_cond = max(cond.values()) if len(cond) else 0

        rows.append(
            {
                "shell_sample_id": r["shell_sample_id"],
                "k": r["k"],
                "domain_neighbor_purity": (max_dom / n_dom) if n_dom else np.nan,
                "condition_neighbor_purity": (max_cond / n_cond) if n_cond else np.nan,
                "dominant_neighbor_domain": max(dom, key=dom.get) if len(dom) else np.nan,
                "dominant_neighbor_condition": max(cond, key=cond.get) if len(cond) else np.nan,
                "neighbor_median_distance": r.get("neighbor_median_distance"),
                "neighbor_phi_median": r.get("neighbor_phi_median"),
            }
        )

    out = pd.DataFrame(rows)
    out.to_csv(outdir / "stable_shell_neighbor_purity.csv", index=False)

    print(f"[ok] wrote neighbor purity table: {outdir / 'stable_shell_neighbor_purity.csv'}")


if __name__ == "__main__":
    main()
