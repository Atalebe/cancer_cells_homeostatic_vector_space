from __future__ import annotations

from pathlib import Path
import json
import re
import urllib.request
import urllib.parse
from html.parser import HTMLParser
import yaml


GEO_BASE = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
GEO_SUPPL_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"


def load_config(config_path: str | Path) -> dict:
    with open(config_path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def accession_to_series_stub(accession: str) -> str:
    """
    GEO series FTP layout uses a family folder like GSE124nnn for GSE124989.
    """
    m = re.fullmatch(r"GSE(\d+)", accession.upper())
    if not m:
        raise ValueError(f"Unsupported GEO accession format: {accession}")
    digits = m.group(1)
    if len(digits) < 3:
        raise ValueError(f"Accession digits too short: {accession}")
    prefix = digits[:-3] + "nnn"
    return f"GSE{prefix}"


def build_geo_series_url(accession: str) -> str:
    return f"{GEO_BASE}?acc={urllib.parse.quote(accession)}"


def build_geo_supp_url(accession: str) -> str:
    stub = accession_to_series_stub(accession)
    return f"{GEO_SUPPL_BASE}/{stub}/{accession}/suppl/"


def http_get_text(url: str, timeout: int = 60) -> str:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Mozilla/5.0",
        },
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        charset = resp.headers.get_content_charset() or "utf-8"
        return resp.read().decode(charset, errors="replace")


def http_download(url: str, output_path: str | Path, timeout: int = 120) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Mozilla/5.0",
        },
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp, open(output_path, "wb") as fh:
        fh.write(resp.read())


class HrefCollector(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "a":
            return
        for key, value in attrs:
            if key.lower() == "href" and value:
                self.hrefs.append(value)


def parse_supplementary_links(html: str, accession: str) -> list[str]:
    parser = HrefCollector()
    parser.feed(html)
    out: list[str] = []
    for href in parser.hrefs:
        if accession in href or href.endswith((".txt.gz", ".txt", ".csv", ".tsv", ".mtx.gz", ".mtx", ".gz")):
            if href.startswith("/"):
                href = urllib.parse.urljoin("https://www.ncbi.nlm.nih.gov", href)
            elif href.startswith("ftp://"):
                pass
            elif href.startswith("http://") or href.startswith("https://"):
                pass
            else:
                href = urllib.parse.urljoin(build_geo_supp_url(accession), href)
            out.append(href)
    # preserve order, remove duplicates
    seen = set()
    deduped = []
    for x in out:
        if x not in seen:
            seen.add(x)
            deduped.append(x)
    return deduped


def write_json(path: str | Path, payload: dict) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)


def dry_run_download_plan(config_path: str | Path) -> dict:
    cfg = load_config(config_path)

    raw_dir = ensure_dir(cfg["paths"]["raw_dir"])
    processed_dir = ensure_dir(cfg["paths"]["processed_dir"])
    metadata_dir = ensure_dir(cfg["paths"]["metadata_dir"])

    accession = cfg["dataset"]["accession"]
    series_url = build_geo_series_url(accession)
    supp_url = build_geo_supp_url(accession)

    plan = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "accession": accession,
        "source_type": cfg["dataset"].get("source_type"),
        "title": cfg["dataset"].get("title"),
        "raw_dir": str(raw_dir),
        "processed_dir": str(processed_dir),
        "metadata_dir": str(metadata_dir),
        "download_strategy": cfg["download"].get("strategy"),
        "expected_files": cfg["download"].get("expected_files", []),
        "series_url": series_url,
        "supplementary_url": supp_url,
        "notes": cfg["download"].get("notes", ""),
    }

    write_json(raw_dir / "download_plan.json", plan)
    return plan


def execute_geo_series_download(config_path: str | Path, download_files: bool = False) -> dict:
    cfg = load_config(config_path)
    accession = cfg["dataset"]["accession"]
    raw_dir = ensure_dir(cfg["paths"]["raw_dir"])
    metadata_dir = ensure_dir(cfg["paths"]["metadata_dir"])

    series_url = build_geo_series_url(accession)
    supp_url = build_geo_supp_url(accession)

    series_html = http_get_text(series_url)
    supp_html = http_get_text(supp_url)

    series_html_path = raw_dir / "geo_series_page.html"
    supp_html_path = raw_dir / "geo_supplementary_index.html"

    series_html_path.write_text(series_html, encoding="utf-8")
    supp_html_path.write_text(supp_html, encoding="utf-8")

    supplementary_links = parse_supplementary_links(supp_html, accession)

    manifest = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "accession": accession,
        "series_url": series_url,
        "supplementary_url": supp_url,
        "series_html": str(series_html_path),
        "supplementary_index_html": str(supp_html_path),
        "supplementary_links": supplementary_links,
        "downloaded_files": [],
    }

    if download_files:
        for url in supplementary_links:
            filename = url.rstrip("/").split("/")[-1]
            if not filename:
                continue
            out_path = raw_dir / filename
            try:
                http_download(url, out_path)
                manifest["downloaded_files"].append(str(out_path))
            except Exception as exc:
                manifest.setdefault("download_errors", []).append(
                    {"url": url, "error": str(exc)}
                )

    write_json(raw_dir / "download_manifest.json", manifest)

    metadata_stub = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "accession": accession,
        "expected_populations": cfg.get("metadata", {}).get("expected_populations", []),
        "label_column": cfg.get("metadata", {}).get("label_column"),
        "notes": "Populate after inspecting downloaded supplementary files and GEO page content.",
    }
    write_json(metadata_dir / "metadata_stub.json", metadata_stub)

    return manifest
