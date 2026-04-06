#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary-json", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    summary_path = Path(args.summary_json)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(summary_path, "r", encoding="utf-8") as f:
        payload = json.load(f)

    downloads = payload.get("downloads", [])
    failures = payload.get("failures", [])

    rows = []
    for rec in downloads:
        p = rec.get("path")
        dp = rec.get("decompressed_path")
        row = {
            "url": rec.get("url"),
            "label": rec.get("label"),
            "kind": rec.get("kind"),
            "path": p,
            "path_exists": Path(p).exists() if p else False,
            "path_size_bytes": Path(p).stat().st_size if p and Path(p).exists() else None,
            "decompressed_path": dp,
            "decompressed_exists": Path(dp).exists() if dp else False,
            "decompressed_size_bytes": Path(dp).stat().st_size if dp and Path(dp).exists() else None,
        }
        rows.append(row)

    pd.DataFrame(rows).to_csv(outdir / "bundle_download_records_audit.csv", index=False)
    pd.DataFrame(failures).to_csv(outdir / "bundle_download_failures.csv", index=False)

    print(f"[ok] wrote audit csv: {outdir / 'bundle_download_records_audit.csv'}")
    print(f"[ok] wrote failures csv: {outdir / 'bundle_download_failures.csv'}")
    print(f"[info] n_download_records: {len(rows)}")
    print(f"[info] n_failures: {len(failures)}")


if __name__ == "__main__":
    main()
