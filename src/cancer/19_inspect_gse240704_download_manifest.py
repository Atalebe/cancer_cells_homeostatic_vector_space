from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inspect GSE240704 download manifest and explode useful nested fields."
    )
    parser.add_argument("--manifest-json", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def flatten_links(obj, parent_key=""):
    rows = []

    if isinstance(obj, dict):
        for k, v in obj.items():
            key = f"{parent_key}.{k}" if parent_key else str(k)
            rows.extend(flatten_links(v, key))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            key = f"{parent_key}[{i}]"
            rows.extend(flatten_links(v, key))
    else:
        rows.append({"key": parent_key, "value": obj})

    return rows


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(args.manifest_json, "r", encoding="utf-8", errors="replace") as f:
        manifest = json.load(f)

    with open(args.summary_json, "r", encoding="utf-8", errors="replace") as f:
        summary = json.load(f)

    manifest_flat = pd.DataFrame(flatten_links(manifest))
    summary_flat = pd.DataFrame(flatten_links(summary))

    manifest_flat.to_csv(outdir / "download_manifest_flat.csv", index=False)
    summary_flat.to_csv(outdir / "download_summary_flat.csv", index=False)

    url_like = []
    for df_name, df in [("manifest", manifest_flat), ("summary", summary_flat)]:
        if "value" not in df.columns:
            continue
        sub = df[df["value"].astype(str).str.contains(r"http|ftp|www|geo", case=False, na=False)].copy()
        if len(sub):
            sub.insert(0, "source", df_name)
            url_like.append(sub)

    if url_like:
        pd.concat(url_like, ignore_index=True).to_csv(
            outdir / "url_like_entries.csv",
            index=False,
        )
        print("[ok] wrote url-like entries:", outdir / "url_like_entries.csv")
    else:
        print("[warn] no url-like entries found")

    print("[ok] wrote flat manifest tables to", outdir)
    print("[info] manifest keys:", len(manifest_flat))
    print("[info] summary keys:", len(summary_flat))


if __name__ == "__main__":
    main()
