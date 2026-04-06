#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd
import requests


def ftp_to_https(url: str) -> str:
    if url.startswith("ftp://ftp.ncbi.nlm.nih.gov/"):
        return url.replace("ftp://ftp.ncbi.nlm.nih.gov/", "https://ftp.ncbi.nlm.nih.gov/")
    return url


def safe_name_from_url(url: str, i: int) -> str:
    parsed = urlparse(url)
    tail = Path(parsed.path).name or f"payload_{i}"
    tail = tail.replace("?", "_").replace("&", "_").replace("=", "_")
    if not tail:
        tail = f"payload_{i}"
    return f"{i:02d}_{tail}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidate-csv", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--summary-json", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.candidate_csv)
    records = []

    sess = requests.Session()
    headers = {"User-Agent": "Mozilla/5.0"}

    for i, row in df.iterrows():
        orig_url = str(row["url"])
        url = ftp_to_https(orig_url)
        fname = safe_name_from_url(url, i)
        outpath = outdir / fname

        rec = {
            "orig_url": orig_url,
            "resolved_url": url,
            "outpath": str(outpath),
            "ok": False,
            "status_code": None,
            "content_type": None,
            "size_bytes": None,
            "error": None,
        }

        try:
            r = sess.get(url, headers=headers, timeout=60)
            rec["status_code"] = r.status_code
            rec["content_type"] = r.headers.get("Content-Type")
            if r.ok:
                outpath.write_bytes(r.content)
                rec["ok"] = True
                rec["size_bytes"] = outpath.stat().st_size
            else:
                rec["error"] = f"http_{r.status_code}"
        except Exception as e:
            rec["error"] = repr(e)

        records.append(rec)

    outcsv = outdir / "downloaded_candidate_payloads.csv"
    pd.DataFrame(records).to_csv(outcsv, index=False)

    summary = {
        "candidate_csv": args.candidate_csv,
        "n_attempted": len(records),
        "n_ok": sum(1 for r in records if r["ok"]),
        "n_failed": sum(1 for r in records if not r["ok"]),
        "download_csv": str(outcsv),
    }
    Path(args.summary_json).write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"[ok] wrote payload download table: {outcsv}")
    print(f"[ok] wrote summary: {args.summary_json}")
    print(f"[info] n_attempted: {summary['n_attempted']}")
    print(f"[info] n_ok: {summary['n_ok']}")
    print(f"[info] n_failed: {summary['n_failed']}")


if __name__ == "__main__":
    main()
