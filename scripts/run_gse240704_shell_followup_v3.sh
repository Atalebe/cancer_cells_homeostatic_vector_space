#!/usr/bin/env bash
set -euo pipefail

echo "[followup3] rebuilding selected probe matrix"
PYTHONPATH=. python src/cancer/33_rebuild_gse240704_selected_probe_matrix.py \
  --beta-matrix data/processed/gse240704/beta_matrix.parquet \
  --sample-registry data/metadata/gse240704/sample_registry.csv \
  --outparquet data/processed/gse240704/selected_probe_matrix.parquet \
  --n-probes 5000

echo "[followup3] exporting top probe loadings"
PYTHONPATH=. python src/cancer/31_export_gse240704_top_probe_loadings.py \
  --selected-probe-matrix data/processed/gse240704/selected_probe_matrix.parquet \
  --sample-annotations-curated data/processed/gse240704/sample_annotations_curated.parquet \
  --outdir results/gse240704/top_probe_loadings \
  --top-n 100 \
  --n-components 5

echo "[followup3] testing extreme domain signatures"
PYTHONPATH=. python src/cancer/34_domain_extreme_signature_gse240704.py \
  --merged-curated-csv results/gse240704/domain_condition_followup/domain_condition_merged_curated.csv \
  --outdir results/gse240704/domain_extreme_signature

echo "[followup3] shell neighborhood purity"
PYTHONPATH=. python src/cancer/35_stable_shell_neighbor_purity_gse240704.py \
  --neighbor-summary-csv results/gse240704/stable_shell_neighborhood/stable_shell_neighbor_summary.csv \
  --outdir results/gse240704/stable_shell_neighborhood

echo
echo "[followup3] inspect:"
echo "  results/gse240704/top_probe_loadings/top_probe_loadings_by_pc.csv"
echo "  results/gse240704/top_probe_loadings/pca_explained_variance_rebuilt.csv"
echo "  results/gse240704/domain_extreme_signature/domain_extreme_signature_tests.csv"
echo "  results/gse240704/stable_shell_neighborhood/stable_shell_neighbor_purity.csv"
