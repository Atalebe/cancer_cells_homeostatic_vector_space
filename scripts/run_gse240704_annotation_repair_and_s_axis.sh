#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

SOFT="data/raw/annotations/GSE240704_family.soft.gz"
OUT="results/gse240704/annotation_repair"
mkdir -p "$OUT"

echo "[repair] extracting platform annotation from family soft"
PYTHONPATH=. python src/cancer/45_extract_gpl23976_annotation_from_soft.py \
  --soft-gz "$SOFT" \
  --outcsv "$OUT/GPL23976_annotation_from_soft.csv" \
  --outparquet "$OUT/GPL23976_annotation_from_soft.parquet" \
  --outjson "$OUT/GPL23976_annotation_from_soft_summary.json"

echo "[repair] inspecting parsed platform annotation"
PYTHONPATH=. python src/cancer/46_inspect_parsed_platform_annotation_gse240704.py \
  --annotation-csv "$OUT/GPL23976_annotation_from_soft.csv" \
  --outdir "$OUT/inspection"

echo "[repair] rebuilding PCA sample scores from selected probe matrix"
PYTHONPATH=. python src/cancer/47_rebuild_gse240704_pca_sample_scores.py \
  --selected-probe-matrix data/processed/gse240704/selected_probe_matrix.parquet \
  --outcsv results/gse240704/s_axis_drivers/rebuilt_pca_sample_scores.csv \
  --outjson results/gse240704/s_axis_drivers/rebuilt_pca_sample_scores_summary.json \
  --n-components 5

echo "[repair] mapping S axis to rebuilt PCA and annotating drivers"
PYTHONPATH=. python src/cancer/48_map_s_axis_to_rebuilt_pca_gse240704.py \
  --state-table data/processed/gse240704/state_table.parquet \
  --pca-scores-csv results/gse240704/s_axis_drivers/rebuilt_pca_sample_scores.csv \
  --top-probe-loadings-csv results/gse240704/top_probe_loadings/top_probe_loadings_by_pc.csv \
  --annotation-csv "$OUT/GPL23976_annotation_from_soft.csv" \
  --outdir results/gse240704/s_axis_drivers \
  --top-n 200

echo "[repair] summarizing S-axis driver genes"
PYTHONPATH=. python src/cancer/49_summarize_s_axis_driver_genes_gse240704.py \
  --annotated-csv results/gse240704/s_axis_drivers/s_axis_top_probe_drivers_annotated.csv \
  --outcsv results/gse240704/s_axis_drivers/s_axis_driver_gene_summary.csv \
  --outjson results/gse240704/s_axis_drivers/s_axis_driver_gene_summary.json

echo
echo "[repair] inspect:"
echo "  $OUT/GPL23976_annotation_from_soft_summary.json"
echo "  $OUT/inspection/parsed_annotation_probe_like_columns.csv"
echo "  results/gse240704/s_axis_drivers/rebuilt_pca_sample_scores_summary.json"
echo "  results/gse240704/s_axis_drivers/s_axis_pc_correlation_table.csv"
echo "  results/gse240704/s_axis_drivers/s_axis_driver_summary.json"
echo "  results/gse240704/s_axis_drivers/s_axis_top_probe_drivers_annotated.csv"
echo "  results/gse240704/s_axis_drivers/s_axis_driver_gene_summary.csv"
