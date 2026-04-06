#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

D3CSV="results/gse240704/directional_biology_followup/D3_vs_not_D3_top_200_probes_with_direction.csv"
ANNPAR="results/gse240704/manual_manifest_annotation/gpl23976_annotation_from_manifest.parquet"
SELPROB="data/processed/gse240704/selected_probe_matrix.parquet"
CENTROMERE="data/reference/cytoband/hg19_centromere_breakpoints.csv"

OUTCHR="results/gse240704/chromosome_enrichment_vs_selected_background"
OUTARM="results/gse240704/chromosome_arm_enrichment_followup"
OUTPLOT="results/gse240704/chromosome_enrichment_vs_selected_background/plots"

mkdir -p logs
mkdir -p "$OUTCHR" "$OUTARM" "$OUTPLOT"

echo "[step] chromosome enrichment vs selected-probe background"
PYTHONPATH=. python src/cancer/67_test_d3_chromosome_enrichment_vs_background_gse240704.py \
  --directional-csv "$D3CSV" \
  --background-annotation-parquet "$ANNPAR" \
  --selected-probe-matrix "$SELPROB" \
  --outdir "$OUTCHR" | tee logs/gse240704_chr_vs_background.log

echo
echo "[step] chromosome arm enrichment"
PYTHONPATH=. python src/cancer/68_assign_probe_chromosome_arm_gse240704.py \
  --directional-csv "$D3CSV" \
  --background-annotation-parquet "$ANNPAR" \
  --selected-probe-matrix "$SELPROB" \
  --centromere-csv "$CENTROMERE" \
  --outdir "$OUTARM" | tee logs/gse240704_chr_arm_enrichment.log

echo
echo "[step] plotting chromosome enrichment"
PYTHONPATH=. python src/cancer/69_plot_d3_chromosome_enrichment_gse240704.py \
  --enrichment-csv "$OUTCHR/D3_directional_chromosome_enrichment_vs_selected_background.csv" \
  --outdir "$OUTPLOT" | tee logs/gse240704_chr_enrichment_plot.log

echo
echo "[done] inspect:"
echo "  $OUTCHR/D3_directional_chromosome_enrichment_vs_selected_background.csv"
echo "  $OUTCHR/D3_directional_chromosome_enrichment_vs_selected_background_summary.json"
echo "  $OUTARM/D3_directional_chromosome_arm_enrichment_combined.csv"
echo "  $OUTARM/D3_directional_chromosome_arm_enrichment_summary.json"
echo "  $OUTPLOT/D3_directional_chromosome_log2OR_vs_selected_background.png"
echo "  $OUTPLOT/D3_directional_chromosome_neglog10p_vs_selected_background.png"
echo "  $OUTPLOT/D3_chromosome21_directional_log2OR.png"
