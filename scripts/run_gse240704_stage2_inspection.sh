#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

OUT_BASE="results/gse240704/inspection"
mkdir -p "${OUT_BASE}"

echo "[inspect] stable shell samples"
PYTHONPATH=. python src/cancer/11_inspect_gse240704_stable_shell_samples.py \
  --state-table data/processed/gse240704/state_table.parquet \
  --sample-annotations data/processed/gse240704/sample_annotations.parquet \
  --candidate-shell-table results/gse240704/irreversible_core/candidate_irreversible_shell_table.parquet \
  --state-domains results/gse240704/state_domains/state_domains.parquet \
  --outdir "${OUT_BASE}/stable_shell"

echo "[inspect] state domain centroids"
PYTHONPATH=. python src/cancer/12_inspect_gse240704_state_domain_centroids.py \
  --state-table data/processed/gse240704/state_table.parquet \
  --state-domains results/gse240704/state_domains/state_domains.parquet \
  --sample-annotations data/processed/gse240704/sample_annotations.parquet \
  --outdir "${OUT_BASE}/state_domains"

echo "[inspect] annotation audit"
PYTHONPATH=. python src/cancer/13_audit_gse240704_annotations.py \
  --sample-annotations data/processed/gse240704/sample_annotations.parquet \
  --outdir "${OUT_BASE}/annotations"

echo "[inspect] metadata mapping template"
PYTHONPATH=. python src/cancer/14_build_gse240704_metadata_mapping_template.py \
  --sample-registry data/metadata/gse240704/sample_registry.csv \
  --sample-annotations data/processed/gse240704/sample_annotations.parquet \
  --outcsv "${OUT_BASE}/annotations/gse240704_metadata_mapping_template.csv"

echo "[done] inspection stage finished"
