import argparse
import gzip
import json
import re
import shutil
import urllib.request
from pathlib import Path

import pandas as pd


DEFAULT_GENCODE_GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/"
    "gencode.v49.annotation.gtf.gz"
)


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def download_file(url: str, output_path: Path, timeout: int = 180) -> None:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "Mozilla/5.0"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp, open(output_path, "wb") as fh:
        shutil.copyfileobj(resp, fh)


def parse_gtf_attributes(attr_text: str) -> dict[str, str]:
    out = {}
    for item in attr_text.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        m = re.match(r'(\S+)\s+"([^"]+)"', item)
        if m:
            out[m.group(1)] = m.group(2)
    return out


def build_gene_map_from_gtf(gtf_gz_path: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(gtf_gz_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            feature = parts[2]
            if feature != "gene":
                continue
            attrs = parse_gtf_attributes(parts[8])
            gene_id = attrs.get("gene_id")
            gene_name = attrs.get("gene_name")
            gene_type = attrs.get("gene_type", "")
            if gene_id:
                rows.append(
                    {
                        "ensembl_id": gene_id,
                        "gene_symbol": gene_name if gene_name else gene_id,
                        "gene_type": gene_type,
                    }
                )

    df = pd.DataFrame(rows).drop_duplicates(subset=["ensembl_id"])
    return df.sort_values("ensembl_id").reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--url",
        default=DEFAULT_GENCODE_GTF_URL,
        help="Direct URL to a GENCODE human GTF.gz file",
    )
    parser.add_argument(
        "--output-dir",
        default="data/external/gene_annotations",
        help="Directory to store downloaded annotation files",
    )
    parser.add_argument(
        "--output-csv",
        default="data/external/gene_annotations/gencode_v49_ensembl_to_symbol.csv",
        help="Output mapping CSV path",
    )
    args = parser.parse_args()

    output_dir = ensure_dir(args.output_dir)
    output_csv = Path(args.output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    gtf_gz_path = output_dir / Path(args.url).name

    print(f"[info] downloading: {args.url}")
    download_file(args.url, gtf_gz_path)
    print(f"[ok] downloaded to {gtf_gz_path}")

    gene_map = build_gene_map_from_gtf(gtf_gz_path)
    gene_map.to_csv(output_csv, index=False)

    summary = {
        "source_url": args.url,
        "downloaded_file": str(gtf_gz_path),
        "output_csv": str(output_csv),
        "n_rows": int(gene_map.shape[0]),
        "n_unique_gene_symbols": int(gene_map["gene_symbol"].nunique()),
        "columns": list(gene_map.columns),
    }

    summary_path = output_csv.with_suffix(".summary.json")
    with open(summary_path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
