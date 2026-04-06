#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

echo "[d2d3] probe-level signatures for D2 and D3"
PYTHONPATH=. python src/cancer/36_domain_probe_signature_gse240704.py \
  --selected-probe-matrix data/processed/gse240704/selected_probe_matrix.parquet \
  --merged-curated-csv results/gse240704/domain_condition_followup/domain_condition_merged_curated.csv \
  --outdir results/gse240704/d2_d3_probe_signatures \
  --top-n 200

echo "[d2d3] D3 state audit"
PYTHONPATH=. python src/cancer/37_d3_state_audit_gse240704.py \
  --selected-probe-matrix data/processed/gse240704/selected_probe_matrix.parquet \
  --merged-curated-csv results/gse240704/domain_condition_followup/domain_condition_merged_curated.csv \
  --stable-shell-csv results/gse240704/inspection/stable_shell/stable_shell_samples_inspection.csv \
  --neighbor-csv results/gse240704/stable_shell_neighborhood/stable_shell_k_nearest_neighbors.csv \
  --neighbor-summary-csv results/gse240704/stable_shell_neighborhood/stable_shell_neighbor_summary.csv \
  --outdir results/gse240704/d3_state_audit

echo "[d2d3] plotting D3 diagnostics"
PYTHONPATH=. python src/cancer/38_plot_d3_diagnostics_gse240704.py \
  --sample-metrics-csv results/gse240704/d3_state_audit/d3_sample_level_metrics.csv \
  --outdir results/gse240704/d3_state_audit/plots

echo
echo "[d2d3] inspect:"
echo "  results/gse240704/d2_d3_probe_signatures/domain_probe_signature_summary.json"
echo "  results/gse240704/d2_d3_probe_signatures/D3_vs_not_D3_top_200_probes.csv"
echo "  results/gse240704/d2_d3_probe_signatures/D3_vs_D2_top_200_probes.csv"
echo "  results/gse240704/d2_d3_probe_signatures/D2_vs_D1_top_200_probes.csv"
echo "  results/gse240704/d3_state_audit/d3_state_audit_summary.json"
echo "  results/gse240704/d3_state_audit/d3_metric_tests.csv"
echo "  results/gse240704/d3_state_audit/d3_mgmt_enrichment.csv"
echo "  results/gse240704/d3_state_audit/d3_proximity_to_stable_shell.csv"
echo "  results/gse240704/d3_state_audit/d3_overlap_with_shell_neighborhoods.csv"
echo "  results/gse240704/d3_state_audit/plots/"
