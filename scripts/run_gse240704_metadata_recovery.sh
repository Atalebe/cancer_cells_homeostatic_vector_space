#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

OUTDIR="results/gse240704/metadata_recovery"
RAW_DIR="data/raw/gse240704"

mkdir -p "${OUTDIR}"

echo "[recover] inspecting manifest"
PYTHONPATH=. python src/cancer/19_inspect_gse240704_download_manifest.py \
  --manifest-json data/raw/gse240704/download_manifest.json \
  --summary-json data/raw/gse240704/gse240704_actual_download_summary.json \
  --outdir "${OUTDIR}"

echo "[recover] downloading GEO series metadata"
PYTHONPATH=. python src/cancer/20_download_gse_series_metadata.py \
  --accession GSE240704 \
  --outdir "${RAW_DIR}"

echo "[recover] extracting sample metadata from series matrix"
PYTHONPATH=. python src/cancer/21_extract_gse240704_metadata_from_series_matrix.py \
  --series-matrix data/raw/gse240704/GSE240704_series_matrix.txt \
  --outcsv "${OUTDIR}/gse240704_series_matrix_sample_metadata.csv"

echo "[recover] preparing mapping template with series metadata"
PYTHONPATH=. python src/cancer/22_prepare_gse240704_metadata_mapping_from_series_matrix.py \
  --series-metadata-csv "${OUTDIR}/gse240704_series_matrix_sample_metadata.csv" \
  --mapping-template-csv results/gse240704/inspection/annotations/gse240704_metadata_mapping_template.csv \
  --outcsv "${OUTDIR}/gse240704_metadata_mapping_template_prefilled.csv"

echo
echo "[recover] inspect these files next:"
echo "  ${OUTDIR}/url_like_entries.csv"
echo "  ${OUTDIR}/gse240704_series_matrix_sample_metadata.csv"
echo "  ${OUTDIR}/gse240704_metadata_mapping_template_prefilled.csv"
echo
echo "[recover] if sample_id count and order look sane, copy curated biology fields into:"
echo "  results/gse240704/inspection/annotations/gse240704_metadata_mapping_template.csv"
echo
echo "[recover] then merge with:"
echo "PYTHONPATH=. python src/cancer/18_merge_gse240704_curated_annotations.py \\"
echo "  --sample-annotations data/processed/gse240704/sample_annotations.parquet \\"
echo "  --curated-mapping results/gse240704/inspection/annotations/gse240704_metadata_mapping_template.csv \\"
echo "  --outparquet data/processed/gse240704/sample_annotations_curated.parquet \\"
echo "  --outcsv data/processed/gse240704/sample_annotations_curated.csv"
