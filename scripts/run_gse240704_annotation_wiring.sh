#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

ANN="data/raw/annotations/GPL23976_annotation.parquet"
OUT="results/gse240704/annotation_wiring"
mkdir -p "$OUT"

echo "[annot] inspecting annotation reference"
PYTHONPATH=. python src/cancer/42_inspect_annotation_reference_gse240704.py \
  --annotation "$ANN" \
  --outdir "$OUT/inspection"

echo
echo "[annot] inspect these first:"
echo "  $OUT/inspection/annotation_columns.csv"
echo "  $OUT/inspection/annotation_probe_like_columns.csv"
echo "  $OUT/inspection/annotation_preview_head30.csv"
echo "  $OUT/inspection/annotation_text_columns_summary.csv"

echo
echo "[annot] annotating D3 contrast tables"
PYTHONPATH=. python src/cancer/43_annotate_top_probes_from_parquet_gse240704.py \
  --contrast-csv results/gse240704/d2_d3_probe_signatures/D3_vs_not_D3_top_200_probes.csv \
  --annotation-parquet "$ANN" \
  --outcsv results/gse240704/d3_interpretation_followup/D3_vs_not_D3_top_200_probes_annotated_from_parquet.csv \
  --outjson results/gse240704/d3_interpretation_followup/D3_vs_not_D3_annotation_from_parquet_summary.json

PYTHONPATH=. python src/cancer/43_annotate_top_probes_from_parquet_gse240704.py \
  --contrast-csv results/gse240704/d2_d3_probe_signatures/D3_vs_D2_top_200_probes.csv \
  --annotation-parquet "$ANN" \
  --outcsv results/gse240704/d3_interpretation_followup/D3_vs_D2_top_200_probes_annotated_from_parquet.csv \
  --outjson results/gse240704/d3_interpretation_followup/D3_vs_D2_annotation_from_parquet_summary.json

PYTHONPATH=. python src/cancer/43_annotate_top_probes_from_parquet_gse240704.py \
  --contrast-csv results/gse240704/d2_d3_probe_signatures/D2_vs_D1_top_200_probes.csv \
  --annotation-parquet "$ANN" \
  --outcsv results/gse240704/d3_interpretation_followup/D2_vs_D1_top_200_probes_annotated_from_parquet.csv \
  --outjson results/gse240704/d3_interpretation_followup/D2_vs_D1_annotation_from_parquet_summary.json

echo "[annot] mapping S-axis to PCA and exporting annotated drivers"
PYTHONPATH=. python src/cancer/44_map_s_axis_to_pca_and_export_drivers_gse240704.py \
  --state-table data/processed/gse240704/state_table.parquet \
  --top-probe-loadings results/gse240704/top_probe_loadings/top_probe_loadings_by_pc.csv \
  --annotation-parquet "$ANN" \
  --outdir results/gse240704/s_axis_drivers \
  --top-n 200

echo
echo "[annot] inspect:"
echo "  $OUT/inspection/annotation_probe_like_columns.csv"
echo "  results/gse240704/d3_interpretation_followup/D3_vs_not_D3_top_200_probes_annotated_from_parquet.csv"
echo "  results/gse240704/d3_interpretation_followup/D3_vs_D2_top_200_probes_annotated_from_parquet.csv"
echo "  results/gse240704/d3_interpretation_followup/D2_vs_D1_top_200_probes_annotated_from_parquet.csv"
echo "  results/gse240704/s_axis_drivers/s_axis_pc_correlation_table.csv"
echo "  results/gse240704/s_axis_drivers/s_axis_top_probe_drivers_annotated.csv"
echo "  results/gse240704/s_axis_drivers/s_axis_driver_summary.json"
