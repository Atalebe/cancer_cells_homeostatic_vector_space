#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Dict, List
from urllib.parse import urljoin, urlparse
from urllib.request import Request, urlopen

USER_AGENT = "Mozilla/5.0"

def fetch_text(url: str, timeout: int = 60) -> str:
    req = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")

def download_file(url: str, outpath: Path, timeout: int = 120) -> Dict[str, object]:
    req = Request(url, headers={"User-Agent": USER_AGENT})
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(req, timeout=timeout) as resp, open(outpath, "wb") as f:
        shutil.copyfileobj(resp, f)
    return {
        "url": url,
        "path": str(outpath),
        "size_bytes": outpath.stat().st_size,
    }

def sanitize_filename(name: str) -> str:
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    return name or "downloaded_file"

def extract_candidate_links(html: str, base_url: str) -> List[Dict[str, str]]:
    # Very plain parser, enough for GEO accession pages.
    matches = re.findall(r'<a[^>]+href="([^"]+)"[^>]*>(.*?)</a>', html, flags=re.I | re.S)
    out: List[Dict[str, str]] = []
    seen = set()

    for href, label in matches:
        clean_label = re.sub(r"<.*?>", "", label).strip()
        full = urljoin(base_url, href)
        text = f"{clean_label} {full}".lower()

        keep = False
        if "soft" in text:
            keep = True
        if "miniml" in text:
            keep = True
        if "supplementary" in text:
            keep = True
        if "download" in text and "gpl23976" in text:
            keep = True
        if "geo/query/acc.cgi" in full and "gpl23976" in full.lower():
            keep = True

        if keep and full not in seen:
            seen.add(full)
            out.append({"label": clean_label, "url": full})
    return out

def maybe_decompress_gzip(path: Path) -> str | None:
    if path.suffix != ".gz":
        return None
    outpath = path.with_suffix("")
    try:
        with gzip.open(path, "rb") as fin, open(outpath, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        return str(outpath)
    except Exception:
        return None

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--summary-json", required=True)
    args = ap.parse_args()

    accession = args.accession.strip()
    outdir = Path(args.outdir)
    summary_json = Path(args.summary_json)
    outdir.mkdir(parents=True, exist_ok=True)
    summary_json.parent.mkdir(parents=True, exist_ok=True)

    page_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
    html = fetch_text(page_url)

    accession_html = outdir / f"{accession}_accession_page.html"
    accession_html.write_text(html, encoding="utf-8")

    links = extract_candidate_links(html, page_url)

    downloads: List[Dict[str, object]] = []
    failures: List[Dict[str, str]] = []

    # Always keep the accession page.
    downloads.append({
        "url": page_url,
        "path": str(accession_html),
        "size_bytes": accession_html.stat().st_size,
        "kind": "accession_html",
    })

    for item in links:
        url = item["url"]
        parsed = urlparse(url)
        fname = sanitize_filename(Path(parsed.path).name or item["label"] or "downloaded_resource")
        if not fname:
            fname = "downloaded_resource"
        outpath = outdir / fname
        if outpath.exists() and outpath.stat().st_size > 0:
            rec = {
                "url": url,
                "path": str(outpath),
                "size_bytes": outpath.stat().st_size,
                "kind": "cached",
                "label": item["label"],
            }
            maybe_plain = maybe_decompress_gzip(outpath)
            if maybe_plain is not None:
                rec["decompressed_path"] = maybe_plain
            downloads.append(rec)
            continue
        try:
            rec = download_file(url, outpath)
            rec["kind"] = "downloaded"
            rec["label"] = item["label"]
            maybe_plain = maybe_decompress_gzip(outpath)
            if maybe_plain is not None:
                rec["decompressed_path"] = maybe_plain
            downloads.append(rec)
        except Exception as e:
            failures.append({"url": url, "label": item["label"], "error": repr(e)})

    summary = {
        "accession": accession,
        "accession_page": page_url,
        "outdir": str(outdir),
        "n_candidate_links": len(links),
        "n_download_records": len(downloads),
        "n_failures": len(failures),
        "downloads": downloads,
        "failures": failures,
    }
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"[ok] wrote download summary: {summary_json}")
    print(f"[info] accession: {accession}")
    print(f"[info] candidate links: {len(links)}")
    print(f"[info] download records: {len(downloads)}")
    print(f"[info] failures: {len(failures)}")

if __name__ == "__main__":
    main()
