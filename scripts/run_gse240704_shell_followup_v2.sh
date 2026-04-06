#!/usr/bin/env bash
set -euo pipefail

MERGED="results/gse240704/domain_condition_followup/domain_condition_merged_curated.csv"

if [[ ! -f "$MERGED" ]]; then
  echo "[followup2] merged curated table missing, rebuilding it"
  PYTHONPATH=. python src/cancer/32_build_gse240704_domain_condition_merged_curated.py \
    --state-table data/processed/gse240704/state_table.parquet \
    --sample-annotations-curated data/processed/gse240704/sample_annotations_curated.parquet \
    --state-domains-csv results/gse240704/state_domains/state_domains.csv \
    --stable-shell-csv results/gse240704/inspection/stable_shell/stable_shell_samples_inspection.csv \
    --outdir results/gse240704/domain_condition_followup
fi

echo "[followup2] pairwise domain contrasts"
PYTHONPATH=. python src/cancer/29_pairwise_domain_contrasts_gse240704.py \
  --domain-followup-dir results/gse240704/domain_condition_followup \
  --outdir results/gse240704/domain_pairwise_followup

echo "[followup2] stable-shell neighborhood analysis"
PYTHONPATH=. python src/cancer/30_stable_shell_neighborhood_gse240704.py \
  --merged-curated-csv results/gse240704/domain_condition_followup/domain_condition_merged_curated.csv \
  --stable-shell-csv results/gse240704/inspection/stable_shell/stable_shell_samples_inspection.csv \
  --outdir results/gse240704/stable_shell_neighborhood \
  --k 15

echo "[followup2] top probe loading export"
SELECTED_MATRIX=""
for CAND in \
  data/processed/gse240704/selected_probe_matrix.parquet \
  data/processed/gse240704/selected_probe_matrix.csv \
  data/processed/gse240704/hrsm_selected_probe_matrix.parquet \
  data/processed/gse240704/hrsm_selected_probe_matrix.csv \
  results/gse240704/hrsm_selected_probe_matrix.parquet \
  results/gse240704/hrsm_selected_probe_matrix.csv
do
  if [[ -f "$CAND" ]]; then
    SELECTED_MATRIX="$CAND"
    break
  fi
done

if [[ -n "$SELECTED_MATRIX" ]]; then
  PYTHONPATH=. python src/cancer/31_export_gse240704_top_probe_loadings.py \
    --selected-probe-matrix "$SELECTED_MATRIX" \
    --sample-annotations-curated data/processed/gse240704/sample_annotations_curated.parquet \
    --outdir results/gse240704/top_probe_loadings \
    --top-n 100 \
    --n-components 5
else
  echo "[warn] no saved selected probe matrix found, skipping top-probe loading export"
fi

echo
echo "[followup2] inspect:"
echo "  results/gse240704/domain_condition_followup/domain_condition_merged_curated.csv"
echo "  results/gse240704/domain_pairwise_followup/pairwise_domain_contrasts.csv"
echo "  results/gse240704/stable_shell_neighborhood/stable_shell_neighbor_summary.csv"
echo "  results/gse240704/stable_shell_neighborhood/stable_shell_k_nearest_neighbors.csv"
echo "  results/gse240704/top_probe_loadings/top_probe_loadings_by_pc.csv"
