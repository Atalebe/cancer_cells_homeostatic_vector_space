#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

echo "[merge] building MGMT curated mapping"
PYTHONPATH=. python src/cancer/23_build_gse240704_mgmt_mapping.py \
  --series-metadata-csv results/gse240704/metadata_recovery/gse240704_series_matrix_sample_metadata.csv \
  --mapping-template-csv results/gse240704/inspection/annotations/gse240704_metadata_mapping_template.csv \
  --outcsv results/gse240704/inspection/annotations/gse240704_metadata_mapping_template_curated_auto.csv

echo "[merge] promoting curated mapping into main inspection path"
cp results/gse240704/inspection/annotations/gse240704_metadata_mapping_template_curated_auto.csv \
   results/gse240704/inspection/annotations/gse240704_metadata_mapping_template.csv

echo "[merge] merging curated annotations"
PYTHONPATH=. python src/cancer/18_merge_gse240704_curated_annotations.py \
  --sample-annotations data/processed/gse240704/sample_annotations.parquet \
  --curated-mapping results/gse240704/inspection/annotations/gse240704_metadata_mapping_template.csv \
  --outparquet data/processed/gse240704/sample_annotations_curated.parquet \
  --outcsv data/processed/gse240704/sample_annotations_curated.csv

echo "[merge] rerunning curated condition-axis tests"
PYTHONPATH=. python src/cancer/24_rerun_gse240704_condition_axis_tests_curated.py \
  --state-table data/processed/gse240704/state_table.parquet \
  --sample-annotations-curated data/processed/gse240704/sample_annotations_curated.parquet \
  --stable-shell-csv results/gse240704/inspection/stable_shell/stable_shell_samples_inspection.csv \
  --outdir results/gse240704/condition_axis_tests_curated

echo
echo "[done] inspect:"
echo "  results/gse240704/inspection/annotations/gse240704_metadata_mapping_template.csv"
echo "  data/processed/gse240704/sample_annotations_curated.csv"
echo "  results/gse240704/condition_axis_tests_curated/condition_axis_tests_curated.csv"
echo "  results/gse240704/condition_axis_tests_curated/stable_shell_mgmt_overlay.csv"
echo "  results/gse240704/condition_axis_tests_curated/condition_axis_summary.json"
