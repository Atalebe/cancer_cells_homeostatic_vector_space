#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

RAW_DIR="data/raw/gse240704"
OUTDIR="results/gse240704/metadata_inspection"

mkdir -p "${OUTDIR}"

echo "[meta] inspecting raw metadata sources"
PYTHONPATH=. python src/cancer/15_inspect_gse240704_raw_metadata_sources.py \
  --raw-dir "${RAW_DIR}" \
  --outdir "${OUTDIR}"

echo
echo "[meta] done. now inspect:"
echo "  ${OUTDIR}/raw_metadata_source_inventory.csv"
echo "  ${OUTDIR}/raw_metadata_candidate_columns_union.csv"
echo
echo "[meta] choose the best source file, then run:"
echo "PYTHONPATH=. python src/cancer/16_extract_gse240704_sample_metadata_table.py \\"
echo "  --input <chosen_metadata_file> \\"
echo "  --outcsv ${OUTDIR}/selected_sample_metadata_table.csv"
echo
echo "PYTHONPATH=. python src/cancer/17_match_gse240704_metadata_to_sample_ids.py \\"
echo "  --sample-registry data/metadata/gse240704/sample_registry.csv \\"
echo "  --metadata-csv ${OUTDIR}/selected_sample_metadata_table.csv \\"
echo "  --outcsv ${OUTDIR}/matched_sample_metadata_table.csv"
echo
echo "[meta] then copy matched biology into:"
echo "  results/gse240704/inspection/annotations/gse240704_metadata_mapping_template.csv"
echo
echo "[meta] after manual curation, merge it with:"
echo "PYTHONPATH=. python src/cancer/18_merge_gse240704_curated_annotations.py \\"
echo "  --sample-annotations data/processed/gse240704/sample_annotations.parquet \\"
echo "  --curated-mapping results/gse240704/inspection/annotations/gse240704_metadata_mapping_template.csv \\"
echo "  --outparquet data/processed/gse240704/sample_annotations_curated.parquet \\"
echo "  --outcsv data/processed/gse240704/sample_annotations_curated.csv"
