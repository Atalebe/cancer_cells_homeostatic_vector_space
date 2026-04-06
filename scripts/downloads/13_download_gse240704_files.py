import argparse
import json
from pathlib import Path
from urllib.request import urlretrieve

import yaml


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--include-raw-tar", action="store_true")
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    raw_dir = Path(cfg["paths"]["raw_dir"])
    raw_dir.mkdir(parents=True, exist_ok=True)

    files = {
        "GSE240704_Matrix_non-normalized.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240704/suppl/GSE240704_Matrix_non-normalized.txt.gz",
        "GSE240704_Matrix_normalized.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240704/suppl/GSE240704_Matrix_normalized.txt.gz",
        "filelist.txt": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240704/suppl/filelist.txt",
    }

    if args.include_raw_tar:
        files["GSE240704_RAW.tar"] = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240704/suppl/GSE240704_RAW.tar"

    downloaded = []
    skipped = []

    for fname, url in files.items():
        out = raw_dir / fname
        if out.exists() and out.stat().st_size > 0:
            skipped.append(str(out))
            continue
        print(f"[info] downloading {fname}")
        urlretrieve(url, out)
        downloaded.append(str(out))

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "downloaded_files": downloaded,
        "skipped_existing_files": skipped,
        "raw_dir": str(raw_dir),
    }

    with open(raw_dir / "gse240704_actual_download_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
