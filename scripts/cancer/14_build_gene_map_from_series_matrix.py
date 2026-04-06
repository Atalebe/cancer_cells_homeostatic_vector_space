import argparse
import json
from pathlib import Path

import pandas as pd
import yaml


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--mapping-file", default=None)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    processed_dir.mkdir(parents=True, exist_ok=True)

    expr_path = processed_dir / "expression_matrix_wide.parquet"
    if not expr_path.exists():
        raise FileNotFoundError(f"Missing expression matrix: {expr_path}")

    expr = pd.read_parquet(expr_path)
    gene_col = expr.columns[0]

    out_path = processed_dir / "gene_map.csv"

    if args.mapping_file:
        mapping_file = Path(args.mapping_file)
        if not mapping_file.exists():
            raise FileNotFoundError(f"Missing mapping file: {mapping_file}")
        gm = pd.read_csv(mapping_file)
        gm.to_csv(out_path, index=False)

        summary = {
            "dataset_id": cfg["dataset"]["dataset_id"],
            "source": str(mapping_file),
            "n_rows": int(gm.shape[0]),
            "output_file": str(out_path),
        }
    else:
        fallback = pd.DataFrame(
            {
                "ensembl_id": expr[gene_col].astype(str),
                "gene_symbol": expr[gene_col].astype(str).str.replace(r"\.\d+$", "", regex=True),
            }
        ).drop_duplicates()
        fallback.to_csv(out_path, index=False)

        summary = {
            "dataset_id": cfg["dataset"]["dataset_id"],
            "source": "fallback_from_ensembl_ids",
            "n_rows": int(fallback.shape[0]),
            "output_file": str(out_path),
            "note": "Gene symbols are fallback stripped Ensembl IDs unless a real mapping file is supplied.",
        }

    with open(processed_dir / "gene_map_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
