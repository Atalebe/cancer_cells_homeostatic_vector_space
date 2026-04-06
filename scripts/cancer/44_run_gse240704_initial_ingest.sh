#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

CFG="configs/cancer/gse240704_stem_branch.yaml"

echo "============================================================"
echo "[0/5] download actual supplementary files"
PYTHONPATH=. python scripts/downloads/13_download_gse240704_files.py \
  --config "$CFG"

echo "============================================================"
echo "[1/5] inspect downloads"
PYTHONPATH=. python scripts/cancer/40_inspect_gse240704_downloads.py \
  --config "$CFG"

echo "============================================================"
echo "[2/5] ingest normalized matrix"
PYTHONPATH=. python scripts/cancer/41_ingest_gse240704_matrices.py \
  --config "$CFG" \
  --matrix-kind normalized

echo "============================================================"
echo "[3/5] build metadata registry"
PYTHONPATH=. python scripts/cancer/42_build_gse240704_metadata_registry.py \
  --config "$CFG" \
  --matrix-kind normalized

echo "============================================================"
echo "[4/5] build master table"
PYTHONPATH=. python scripts/cancer/43_build_gse240704_master_table.py \
  --config "$CFG" \
  --matrix-kind normalized

echo "============================================================"
echo "[done] gse240704 initial inspect and ingest complete"
