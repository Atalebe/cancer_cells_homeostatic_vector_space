#!/usr/bin/env python3

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path
from urllib.request import urlretrieve

import yaml


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("[run]", " ".join(cmd))
    subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True)


def main() -> None:
    cfg = read_config()
    raw_dir = REPO_ROOT / cfg["raw_dir"]
    mkdir(raw_dir)

    series = cfg["geo"]["series_accession"]
    superseries = cfg["geo"]["superseries_accession"]

    # GEO endpoints
    soft_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series[:-3]}nnn/{series}/soft/{series}_family.soft.gz"
    miniml_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series[:-3]}nnn/{series}/miniml/{series}_family.xml.tgz"
    matrix_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series[:-3]}nnn/{series}/matrix/{series}_series_matrix.txt.gz"

    superseries_soft_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{superseries[:-3]}nnn/{superseries}/soft/{superseries}_family.soft.gz"
    superseries_matrix_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{superseries[:-3]}nnn/{superseries}/matrix/{superseries}_series_matrix.txt.gz"

    files = {
        "gse161895_family_soft_gz": raw_dir / f"{series}_family.soft.gz",
        "gse161895_family_miniml_tgz": raw_dir / f"{series}_family.xml.tgz",
        "gse161895_series_matrix_gz": raw_dir / f"{series}_series_matrix.txt.gz",
        "gse161901_family_soft_gz": raw_dir / f"{superseries}_family.soft.gz",
        "gse161901_series_matrix_gz": raw_dir / f"{superseries}_series_matrix.txt.gz",
    }

    url_map = {
        "gse161895_family_soft_gz": soft_url,
        "gse161895_family_miniml_tgz": miniml_url,
        "gse161895_series_matrix_gz": matrix_url,
        "gse161901_family_soft_gz": superseries_soft_url,
        "gse161901_series_matrix_gz": superseries_matrix_url,
    }

    manifest = {}
    for key, outpath in files.items():
        url = url_map[key]
        try:
            print(f"[download] {url}")
            urlretrieve(url, outpath)
            manifest[key] = {"url": url, "status": "downloaded", "path": str(outpath)}
        except Exception as e:
            manifest[key] = {"url": url, "status": f"failed: {e}", "path": str(outpath)}

    # Try to grab GEO accession page HTML as a frozen reference
    try:
        html_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={series}"
        html_path = raw_dir / f"{series}_geo_page.html"
        urlretrieve(html_url, html_path)
        manifest["geo_html"] = {"url": html_url, "status": "downloaded", "path": str(html_path)}
    except Exception as e:
        manifest["geo_html"] = {"url": html_url, "status": f"failed: {e}"}

    manifest_path = raw_dir / "download_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    print(f"[ok] wrote {manifest_path}")


if __name__ == "__main__":
    main()
