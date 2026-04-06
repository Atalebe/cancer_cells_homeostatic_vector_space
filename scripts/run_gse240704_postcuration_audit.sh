#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

echo "[post] repairing stable-shell overlay"
PYTHONPATH=. python src/cancer/25_fix_gse240704_stable_shell_overlay.py \
  --overlay-csv results/gse240704/condition_axis_tests_curated/stable_shell_mgmt_overlay.csv \
  --stable-shell-csv results/gse240704/inspection/stable_shell/stable_shell_samples_inspection.csv \
  --sample-annotations-curated data/processed/gse240704/sample_annotations_curated.parquet \
  --outdir results/gse240704/condition_axis_tests_curated

echo "[post] auditing target stable-shell mappings"
PYTHONPATH=. python src/cancer/26_audit_gse240704_sample341_mapping.py

echo "[post] done"
echo "inspect:"
echo "  results/gse240704/condition_axis_tests_curated/stable_shell_mgmt_overlay_repaired.csv"
echo "  results/gse240704/condition_axis_tests_curated/stable_shell_missing_annotation_audit.csv"
echo "  results/gse240704/condition_axis_tests_curated/sample_annotations_curated_target_samples.csv"
