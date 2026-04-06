#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

PYTHONPATH=. python src/cancer/00_download_gse161895.py
PYTHONPATH=. python src/cancer/01_build_gse161895_metadata_registry.py
PYTHONPATH=. python src/cancer/02_audit_gse161895_dataset.py

echo
echo "[next] put the processed expression matrix into data/raw/gse161895/"
echo "[next] accepted names:"
echo "  expression_matrix_wide.parquet"
echo "  expression_matrix_wide.csv"
echo "  expression_matrix_wide.tsv"
echo "  counts_matrix.tsv"
echo "  counts_matrix.csv"
echo "  counts_matrix.parquet"
