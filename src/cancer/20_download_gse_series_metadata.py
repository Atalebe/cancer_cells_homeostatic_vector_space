from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path
from urllib.request import urlopen, Request


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download GEO series matrix and family SOFT files for a GSE accession."
    )
    parser.add_argument("--accession", required=True, help="Example: GSE240704")
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def geo_series_stub(accession: str) -> str:
    digits = accession.replace("GSE", "")
    if len(digits) < 3:
        raise ValueError(f"Unexpected accession: {accession}")
    return f"GSE{digits[:-3]}nnn"


def download(url: str, outpath: Path) -> None:
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req) as r, open(outpath, "wb") as f:
        shutil.copyfileobj(r, f)


def gunzip_file(src: Path, dst: Path) -> None:
    with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout)


def main() -> None:
    args = parse_args()
    accession = args.accession.strip()
    stub = geo_series_stub(accession)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    matrix_url = (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{accession}/matrix/"
        f"{accession}_series_matrix.txt.gz"
    )
    soft_url = (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{accession}/soft/"
        f"{accession}_family.soft.gz"
    )

    matrix_gz = outdir / f"{accession}_series_matrix.txt.gz"
    matrix_txt = outdir / f"{accession}_series_matrix.txt"
    soft_gz = outdir / f"{accession}_family.soft.gz"
    soft_txt = outdir / f"{accession}_family.soft"

    for url, gz_path, txt_path, label in [
        (matrix_url, matrix_gz, matrix_txt, "series matrix"),
        (soft_url, soft_gz, soft_txt, "family soft"),
    ]:
        try:
            print(f"[download] {label}: {url}")
            download(url, gz_path)
            print(f"[ok] downloaded: {gz_path}")
            gunzip_file(gz_path, txt_path)
            print(f"[ok] unpacked: {txt_path}")
        except Exception as e:
            print(f"[warn] failed {label}: {e}")

    print("[done] metadata download attempt finished")


if __name__ == "__main__":
    main()
