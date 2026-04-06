#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

mkdir -p data/raw/annotations/gpl23976_bundle
mkdir -p results/gse240704/annotation_bundle_download

PYTHONPATH=. python src/cancer/50_download_gpl23976_reference_bundle.py \
  --accession GPL23976 \
  --outdir data/raw/annotations/gpl23976_bundle \
  --summary-json results/gse240704/annotation_bundle_download/gpl23976_bundle_download_summary.json

echo
echo "[download] inspect:"
echo "  results/gse240704/annotation_bundle_download/gpl23976_bundle_download_summary.json"
echo "  data/raw/annotations/gpl23976_bundle"
