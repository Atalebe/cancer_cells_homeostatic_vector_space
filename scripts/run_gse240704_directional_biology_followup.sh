#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

OUTDIR="results/gse240704/directional_biology_followup"
mkdir -p "$OUTDIR"

echo "[dir] step 1, add direction to contrast tables"
PYTHONPATH=. python src/cancer/64_add_direction_to_reannotated_tables_gse240704.py \
  --inputs \
    results/gse240704/manifest_reannotated/D3_vs_not_D3_top_200_probes_normalized_annotation.csv \
    results/gse240704/manifest_reannotated/D3_vs_D2_top_200_probes_normalized_annotation.csv \
    results/gse240704/manifest_reannotated/D2_vs_D1_top_200_probes_normalized_annotation.csv \
  --outdir "$OUTDIR"

echo "[dir] step 2, prepare signed S-axis table"
PYTHONPATH=. python src/cancer/67_prepare_s_axis_signed_table_gse240704.py \
  --input results/gse240704/manifest_reannotated/s_axis_top_probe_drivers_annotated_normalized_annotation.csv \
  --outcsv "$OUTDIR/s_axis_top_probe_drivers_with_direction.csv" \
  --summary-json "$OUTDIR/s_axis_top_probe_drivers_with_direction_summary.json"

echo "[dir] step 3, directional gene summaries"
PYTHONPATH=. python src/cancer/65_summarize_directional_gene_tokens_gse240704.py \
  --input "$OUTDIR/D3_vs_not_D3_top_200_probes_with_direction.csv" \
  --label D3_vs_not_D3 \
  --outdir "$OUTDIR"

PYTHONPATH=. python src/cancer/65_summarize_directional_gene_tokens_gse240704.py \
  --input "$OUTDIR/D3_vs_D2_top_200_probes_with_direction.csv" \
  --label D3_vs_D2 \
  --outdir "$OUTDIR"

PYTHONPATH=. python src/cancer/65_summarize_directional_gene_tokens_gse240704.py \
  --input "$OUTDIR/D2_vs_D1_top_200_probes_with_direction.csv" \
  --label D2_vs_D1 \
  --outdir "$OUTDIR"

PYTHONPATH=. python src/cancer/65_summarize_directional_gene_tokens_gse240704.py \
  --input "$OUTDIR/s_axis_top_probe_drivers_with_direction.csv" \
  --label S_axis \
  --outdir "$OUTDIR"

echo "[dir] step 4, directional context summaries"
PYTHONPATH=. python src/cancer/66_summarize_directional_context_features_gse240704.py \
  --input "$OUTDIR/D3_vs_not_D3_top_200_probes_with_direction.csv" \
  --label D3_vs_not_D3 \
  --outdir "$OUTDIR"

PYTHONPATH=. python src/cancer/66_summarize_directional_context_features_gse240704.py \
  --input "$OUTDIR/D3_vs_D2_top_200_probes_with_direction.csv" \
  --label D3_vs_D2 \
  --outdir "$OUTDIR"

PYTHONPATH=. python src/cancer/66_summarize_directional_context_features_gse240704.py \
  --input "$OUTDIR/D2_vs_D1_top_200_probes_with_direction.csv" \
  --label D2_vs_D1 \
  --outdir "$OUTDIR"

PYTHONPATH=. python src/cancer/66_summarize_directional_context_features_gse240704.py \
  --input "$OUTDIR/s_axis_top_probe_drivers_with_direction.csv" \
  --label S_axis \
  --outdir "$OUTDIR"

echo
echo "[dir] inspect:"
echo "  $OUTDIR/D3_vs_not_D3_directional_gene_token_counts.csv"
echo "  $OUTDIR/D3_vs_not_D3_refgene_group_by_direction.csv"
echo "  $OUTDIR/D3_vs_not_D3_relation_to_cpg_island_by_direction.csv"
echo "  $OUTDIR/S_axis_directional_gene_token_counts.csv"
echo "  $OUTDIR/S_axis_relation_to_cpg_island_by_direction.csv"
