#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

OUTDIR="results/gse240704/biology_interpretation_followup"
mkdir -p "$OUTDIR"

echo "[bio] gene symbol summaries"
PYTHONPATH=. python src/cancer/60_summarize_reannotated_gene_symbols_gse240704.py \
  --input-csv results/gse240704/manifest_reannotated/D3_vs_not_D3_top_200_probes_normalized_annotation.csv \
  --outdir "$OUTDIR" \
  --label D3_vs_not_D3

PYTHONPATH=. python src/cancer/60_summarize_reannotated_gene_symbols_gse240704.py \
  --input-csv results/gse240704/manifest_reannotated/D3_vs_D2_top_200_probes_normalized_annotation.csv \
  --outdir "$OUTDIR" \
  --label D3_vs_D2

PYTHONPATH=. python src/cancer/60_summarize_reannotated_gene_symbols_gse240704.py \
  --input-csv results/gse240704/manifest_reannotated/D2_vs_D1_top_200_probes_normalized_annotation.csv \
  --outdir "$OUTDIR" \
  --label D2_vs_D1

PYTHONPATH=. python src/cancer/60_summarize_reannotated_gene_symbols_gse240704.py \
  --input-csv results/gse240704/manifest_reannotated/s_axis_top_probe_drivers_annotated_normalized_annotation.csv \
  --outdir "$OUTDIR" \
  --label S_axis

echo "[bio] context summaries"
PYTHONPATH=. python src/cancer/61_summarize_reannotated_context_features_gse240704.py \
  --input-csv results/gse240704/manifest_reannotated/D3_vs_not_D3_top_200_probes_normalized_annotation.csv \
  --outdir "$OUTDIR" \
  --label D3_vs_not_D3

PYTHONPATH=. python src/cancer/61_summarize_reannotated_context_features_gse240704.py \
  --input-csv results/gse240704/manifest_reannotated/D3_vs_D2_top_200_probes_normalized_annotation.csv \
  --outdir "$OUTDIR" \
  --label D3_vs_D2

PYTHONPATH=. python src/cancer/61_summarize_reannotated_context_features_gse240704.py \
  --input-csv results/gse240704/manifest_reannotated/D2_vs_D1_top_200_probes_normalized_annotation.csv \
  --outdir "$OUTDIR" \
  --label D2_vs_D1

PYTHONPATH=. python src/cancer/61_summarize_reannotated_context_features_gse240704.py \
  --input-csv results/gse240704/manifest_reannotated/s_axis_top_probe_drivers_annotated_normalized_annotation.csv \
  --outdir "$OUTDIR" \
  --label S_axis

echo "[bio] S-axis biology summary"
PYTHONPATH=. python src/cancer/63_summarize_s_axis_driver_biology_gse240704.py \
  --input-csv results/gse240704/manifest_reannotated/s_axis_top_probe_drivers_annotated_normalized_annotation.csv \
  --outdir "$OUTDIR"

echo
echo "[bio] inspect:"
echo "  $OUTDIR/D3_vs_not_D3_gene_symbol_summary.json"
echo "  $OUTDIR/D3_vs_not_D3_gene_token_counts.csv"
echo "  $OUTDIR/D3_vs_not_D3_relation_to_cpg_island_summary.csv"
echo "  $OUTDIR/S_axis_gene_symbol_summary.json"
echo "  $OUTDIR/s_axis_driver_gene_token_counts.csv"
echo "  $OUTDIR/s_axis_driver_relation_to_cpg_island_summary.csv"
echo "  $OUTDIR/s_axis_driver_biology_summary.json"
