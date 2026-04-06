#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

OUTDIR="results/gse240704/final_directional_reports"
CHRDIR="results/gse240704/chromosome_enrichment_followup"
mkdir -p "$OUTDIR" "$CHRDIR"

echo "[step] building D3 compact directional report"
PYTHONPATH=. python src/cancer/64_build_d3_directional_compact_report_gse240704.py \
  --input-csv results/gse240704/directional_biology_followup/D3_vs_not_D3_top_200_probes_with_direction.csv \
  --outdir "$OUTDIR" \
  --label D3_vs_not_D3 \
  --top-n 30

echo "[step] building S-axis positive loading report"
PYTHONPATH=. python src/cancer/65_build_s_axis_positive_loading_report_gse240704.py \
  --input-csv results/gse240704/directional_biology_followup/s_axis_top_probe_drivers_with_direction.csv \
  --outdir "$OUTDIR" \
  --positive-label positive_loading \
  --top-n 40

echo "[step] testing D3 chromosome enrichment"
PYTHONPATH=. python src/cancer/66_test_d3_chromosomal_enrichment_gse240704.py \
  --input-csv results/gse240704/directional_biology_followup/D3_vs_not_D3_top_200_probes_with_direction.csv \
  --outdir "$CHRDIR" \
  --higher-label higher_in_a \
  --lower-label lower_in_a

echo
echo "[done] inspect:"
echo "  $OUTDIR/D3_vs_not_D3_directional_compact_report.csv"
echo "  $OUTDIR/s_axis_positive_loading_top_genes.csv"
echo "  $OUTDIR/s_axis_positive_loading_top_probes.csv"
echo "  $CHRDIR/D3_higher_in_a_chromosome_enrichment.csv"
echo "  $CHRDIR/D3_lower_in_a_chromosome_enrichment.csv"
echo "  $CHRDIR/D3_directional_chromosome_counts.csv"
