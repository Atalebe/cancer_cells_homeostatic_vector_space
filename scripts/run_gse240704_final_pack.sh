#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

echo "[final] exporting compact final result bundle"
PYTHONPATH=. python src/cancer/67_export_gse240704_final_result_bundle.py \
  --base results/gse240704 \
  --outdir results/gse240704/final_result_bundle

echo "[final] building manuscript-ready tables"
PYTHONPATH=. python src/cancer/68_build_gse240704_manuscript_tables.py \
  --base results/gse240704 \
  --outdir results/gse240704/manuscript_tables

echo "[final] plotting integrated summary panels"
PYTHONPATH=. python src/cancer/69_plot_gse240704_final_integrated_panels.py \
  --base results/gse240704 \
  --outdir results/gse240704/final_summary_plots

echo "[final] building candidate headline table"
PYTHONPATH=. python src/cancer/70_build_gse240704_candidate_headline_table.py \
  --base results/gse240704 \
  --outcsv results/gse240704/final_result_bundle/gse240704_candidate_headline_table.csv

echo "[final] writing handover summary"
PYTHONPATH=. python src/cancer/71_build_gse240704_handover_summary.py \
  --base results/gse240704 \
  --outjson results/gse240704/final_result_bundle/gse240704_handover_summary.json

echo
echo "[done] inspect:"
echo "  results/gse240704/final_result_bundle/final_bundle_inventory.csv"
echo "  results/gse240704/final_result_bundle/final_bundle_summary.json"
echo "  results/gse240704/final_result_bundle/gse240704_candidate_headline_table.csv"
echo "  results/gse240704/final_result_bundle/gse240704_handover_summary.json"
echo "  results/gse240704/manuscript_tables/"
echo "  results/gse240704/final_summary_plots/"
