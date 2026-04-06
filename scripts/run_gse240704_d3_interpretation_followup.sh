#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

OUTDIR="results/gse240704/d3_interpretation_followup"
mkdir -p "$OUTDIR"

D3_VS_NOT="results/gse240704/d2_d3_probe_signatures/D3_vs_not_D3_top_200_probes.csv"
D3_VS_D2="results/gse240704/d2_d3_probe_signatures/D3_vs_D2_top_200_probes.csv"
D2_VS_D1="results/gse240704/d2_d3_probe_signatures/D2_vs_D1_top_200_probes.csv"

D3_SAMPLE_LEVEL="results/gse240704/d3_state_audit/d3_sample_level_metrics.csv"
D3_PROXIMITY="results/gse240704/d3_state_audit/d3_proximity_to_stable_shell.csv"

echo "[followup4] annotating top contrast probes"
PYTHONPATH=. python src/cancer/39_annotate_top_probes_gse240704.py \
  --contrast-csv "$D3_VS_NOT" \
  --outcsv "$OUTDIR/D3_vs_not_D3_top_200_probes_annotated.csv" \
  --outjson "$OUTDIR/D3_vs_not_D3_annotation_summary.json"

PYTHONPATH=. python src/cancer/39_annotate_top_probes_gse240704.py \
  --contrast-csv "$D3_VS_D2" \
  --outcsv "$OUTDIR/D3_vs_D2_top_200_probes_annotated.csv" \
  --outjson "$OUTDIR/D3_vs_D2_annotation_summary.json"

PYTHONPATH=. python src/cancer/39_annotate_top_probes_gse240704.py \
  --contrast-csv "$D2_VS_D1" \
  --outcsv "$OUTDIR/D2_vs_D1_top_200_probes_annotated.csv" \
  --outjson "$OUTDIR/D2_vs_D1_annotation_summary.json"

echo "[followup4] summarizing direction-aware signatures"
PYTHONPATH=. python src/cancer/40_summarize_domain_probe_signature_gse240704.py \
  --annotated-contrast-csv "$OUTDIR/D3_vs_not_D3_top_200_probes_annotated.csv" \
  --outcsv "$OUTDIR/D3_vs_not_D3_signature_direction_summary.csv" \
  --outjson "$OUTDIR/D3_vs_not_D3_signature_direction_summary.json" \
  --top-n 200

PYTHONPATH=. python src/cancer/40_summarize_domain_probe_signature_gse240704.py \
  --annotated-contrast-csv "$OUTDIR/D3_vs_D2_top_200_probes_annotated.csv" \
  --outcsv "$OUTDIR/D3_vs_D2_signature_direction_summary.csv" \
  --outjson "$OUTDIR/D3_vs_D2_signature_direction_summary.json" \
  --top-n 200

PYTHONPATH=. python src/cancer/40_summarize_domain_probe_signature_gse240704.py \
  --annotated-contrast-csv "$OUTDIR/D2_vs_D1_top_200_probes_annotated.csv" \
  --outcsv "$OUTDIR/D2_vs_D1_signature_direction_summary.csv" \
  --outjson "$OUTDIR/D2_vs_D1_signature_direction_summary.json" \
  --top-n 200

echo "[followup4] D3 shell-proximity stratification"
PYTHONPATH=. python src/cancer/41_d3_shell_proximity_stratified_audit_gse240704.py \
  --sample-level-csv "$D3_SAMPLE_LEVEL" \
  --proximity-csv "$D3_PROXIMITY" \
  --outcsv "$OUTDIR/d3_shell_proximity_stratified_tests.csv" \
  --outjson "$OUTDIR/d3_shell_proximity_stratified_summary.json"

echo
echo "[followup4] inspect:"
echo "  $OUTDIR/D3_vs_not_D3_annotation_summary.json"
echo "  $OUTDIR/D3_vs_D2_annotation_summary.json"
echo "  $OUTDIR/D2_vs_D1_annotation_summary.json"
echo "  $OUTDIR/D3_vs_not_D3_signature_direction_summary.csv"
echo "  $OUTDIR/D3_vs_D2_signature_direction_summary.csv"
echo "  $OUTDIR/D2_vs_D1_signature_direction_summary.csv"
echo "  $OUTDIR/d3_shell_proximity_stratified_summary.json"
echo "  $OUTDIR/d3_shell_proximity_stratified_tests.csv"
