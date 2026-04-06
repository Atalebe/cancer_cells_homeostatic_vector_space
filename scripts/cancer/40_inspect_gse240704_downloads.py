import argparse
import gzip
import json
from pathlib import Path

import pandas as pd
import yaml


def preview_text_gz(path: Path, n_lines: int = 10) -> list[str]:
    lines = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            lines.append(line.rstrip("\n"))
            if i + 1 >= n_lines:
                break
    return lines


def count_columns_first_data_line(path: Path) -> int:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            return len(line.split("\t"))
    return 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    raw_dir = Path(cfg["paths"]["raw_dir"])
    out_dir = Path(cfg["paths"]["metadata_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    files = {
        "non_normalized": raw_dir / "GSE240704_Matrix_non-normalized.txt.gz",
        "normalized": raw_dir / "GSE240704_Matrix_normalized.txt.gz",
        "raw_tar": raw_dir / "GSE240704_RAW.tar",
        "filelist": raw_dir / "filelist.txt",
    }

    summary = {"dataset_id": cfg["dataset"]["dataset_id"], "files": {}}

    for key, path in files.items():
        entry = {
            "exists": path.exists(),
            "path": str(path),
        }
        if path.exists() and path.suffix == ".gz":
            entry["preview"] = preview_text_gz(path, n_lines=12)
            entry["first_data_line_n_columns"] = count_columns_first_data_line(path)
            entry["size_bytes"] = path.stat().st_size
        elif path.exists():
            entry["size_bytes"] = path.stat().st_size
        summary["files"][key] = entry

    out_path = out_dir / "gse240704_download_inspection.json"
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
