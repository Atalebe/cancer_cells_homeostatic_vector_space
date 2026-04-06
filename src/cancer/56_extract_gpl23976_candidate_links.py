#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from urllib.parse import urljoin

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession-html", required=True)
    ap.add_argument("--outcsv", required=True)
    ap.add_argument("--summary-json", required=True)
    args = ap.parse_args()

    html_path = Path(args.accession_html)
    html = html_path.read_text(encoding="utf-8", errors="replace")

    hrefs = re.findall(r'href=["\']([^"\']+)["\']', html, flags=re.I)
    hrefs = [h.replace("&amp;", "&") for h in hrefs]

    base = "https://www.ncbi.nlm.nih.gov"
    rows = []
    seen = set()

    for href in hrefs:
        full = urljoin(base, href)
        low = full.lower()

        score = 0
        if "gpl23976" in low:
            score += 3
        if any(k in low for k in ["annot", "platform", "table", "txt", "csv", "tsv", "soft", "miniml", "family"]):
            score += 2
        if any(k in low for k in ["ftp.ncbi.nlm.nih.gov", "ncbi.nlm.nih.gov/geo"]):
            score += 1

        if score <= 0:
            continue
        if full in seen:
            continue
        seen.add(full)

        rows.append({
            "url": full,
            "score_hint": score,
            "contains_gpl23976": "gpl23976" in low,
            "contains_annotation_hint": any(k in low for k in ["annot", "platform", "table"]),
            "contains_downloadish_hint": any(k in low for k in ["txt", "csv", "tsv", "soft", "miniml", "family"]),
        })

    df = pd.DataFrame(rows).sort_values(["score_hint", "url"], ascending=[False, True]).reset_index(drop=True)
    Path(args.outcsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.outcsv, index=False)

    summary = {
        "accession_html": str(html_path),
        "n_candidate_links": int(df.shape[0]),
        "top_20": df.head(20).to_dict(orient="records"),
    }
    Path(args.summary_json).write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"[ok] wrote candidate links: {args.outcsv}")
    print(f"[ok] wrote summary: {args.summary_json}")
    print(f"[info] n_candidate_links: {df.shape[0]}")


if __name__ == "__main__":
    main()
