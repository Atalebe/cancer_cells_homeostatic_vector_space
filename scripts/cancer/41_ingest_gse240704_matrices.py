import argparse
import gzip
import json
from pathlib import Path

import pandas as pd
import yaml

from src.common.io import write_table


def read_matrix_gz(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        df = pd.read_csv(fh, sep="\t")
    return df


def normalize_first_column_name(df: pd.DataFrame) -> pd.DataFrame:
    cols = list(df.columns)
    if cols:
        cols[0] = "gene_id"
        df.columns = cols
    return df


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--matrix-kind", choices=["normalized", "non_normalized"], default="normalized")
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    raw_dir = Path(cfg["paths"]["raw_dir"])
    processed_dir = Path(cfg["paths"]["processed_dir"])
    processed_dir.mkdir(parents=True, exist_ok=True)

    file_map = {
        "normalized": raw_dir / "GSE240704_Matrix_normalized.txt.gz",
        "non_normalized": raw_dir / "GSE240704_Matrix_non-normalized.txt.gz",
    }

    matrix_path = file_map[args.matrix_kind]
    if not matrix_path.exists():
        raise FileNotFoundError(f"Missing matrix file: {matrix_path}")

    df = read_matrix_gz(matrix_path)
    df = normalize_first_column_name(df)

    out_parquet = processed_dir / f"gse240704_expression_{args.matrix_kind}.parquet"
    out_csv_head = processed_dir / f"gse240704_expression_{args.matrix_kind}_head.csv"
    write_table(df, out_parquet)
    write_table(df.head(20), out_csv_head)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "matrix_kind": args.matrix_kind,
        "input_file": str(matrix_path),
        "n_rows": int(df.shape[0]),
        "n_columns": int(df.shape[1]),
        "gene_column": df.columns[0],
        "output_parquet": str(out_parquet),
        "output_head_csv": str(out_csv_head),
    }

    with open(processed_dir / f"gse240704_ingest_{args.matrix_kind}_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
