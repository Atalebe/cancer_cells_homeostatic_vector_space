#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

BUNDLE_DIR="data/raw/annotations/gpl23976_bundle"
AUDIT_DIR="results/gse240704/annotation_bundle_download/audit"
INSPECT_DIR="results/gse240704/annotation_bundle_download/inspection"
LINK_DIR="results/gse240704/annotation_bundle_download/link_harvest"
PAYLOAD_DIR="data/raw/annotations/gpl23976_candidate_payloads"
PAYLOAD_INSPECT_DIR="results/gse240704/annotation_bundle_download/payload_inspection"
NORM_DIR="results/gse240704/annotation_bundle_download/normalized"
REANN_DIR="results/gse240704/annotation_bundle_download/reannotated"

echo "[rescue] step 1, download official GPL23976 bundle"
PYTHONPATH=. python src/cancer/50_download_gpl23976_reference_bundle.py \
  --accession GPL23976 \
  --outdir "$BUNDLE_DIR" \
  --summary-json results/gse240704/annotation_bundle_download/gpl23976_bundle_download_summary.json

echo "[rescue] step 2, audit recorded downloads"
PYTHONPATH=. python src/cancer/54_audit_gpl23976_bundle_download_records.py \
  --summary-json results/gse240704/annotation_bundle_download/gpl23976_bundle_download_summary.json \
  --outdir "$AUDIT_DIR"

echo "[rescue] step 3, inspect initial bundle"
PYTHONPATH=. python src/cancer/51_inspect_gpl23976_annotation_bundle.py \
  --indir "$BUNDLE_DIR" \
  --outdir "$INSPECT_DIR"

echo "[rescue] step 4, harvest richer candidate links from accession html"
PYTHONPATH=. python src/cancer/56_extract_gpl23976_candidate_links.py \
  --accession-html "$BUNDLE_DIR/GPL23976_accession_page.html" \
  --outcsv "$LINK_DIR/gpl23976_candidate_links.csv" \
  --summary-json "$LINK_DIR/gpl23976_candidate_links_summary.json"

echo "[rescue] step 5, download candidate payloads"
PYTHONPATH=. python src/cancer/57_download_gpl23976_candidate_payloads.py \
  --candidate-csv "$LINK_DIR/gpl23976_candidate_links.csv" \
  --outdir "$PAYLOAD_DIR" \
  --summary-json "$LINK_DIR/gpl23976_candidate_payload_download_summary.json"

echo "[rescue] step 6, inspect downloaded payloads"
mkdir -p "$PAYLOAD_INSPECT_DIR"
PYTHONPATH=. python src/cancer/51_inspect_gpl23976_annotation_bundle.py \
  --indir "$PAYLOAD_DIR" \
  --outdir "$PAYLOAD_INSPECT_DIR"

echo "[rescue] step 7, auto-pick and normalize best payload if available"
PYTHONPATH=. python src/cancer/55_autopick_and_normalize_gpl23976_annotation.py \
  --inspection-summary "$PAYLOAD_INSPECT_DIR/bundle_inspection_summary.json" \
  --outcsv "$NORM_DIR/gpl23976_annotation_normalized.csv" \
  --outparquet "$NORM_DIR/gpl23976_annotation_normalized.parquet" \
  --summary-json "$NORM_DIR/gpl23976_annotation_normalized_summary.json"

echo "[rescue] step 8, rerun probe annotations"
PYTHONPATH=. python src/cancer/53_rerun_probe_annotations_with_normalized_reference.py \
  --annotation-parquet "$NORM_DIR/gpl23976_annotation_normalized.parquet" \
  --contrast-csvs \
    results/gse240704/d2_d3_probe_signatures/D3_vs_not_D3_top_200_probes.csv \
    results/gse240704/d2_d3_probe_signatures/D3_vs_D2_top_200_probes.csv \
    results/gse240704/d2_d3_probe_signatures/D2_vs_D1_top_200_probes.csv \
  --s-axis-driver-csv results/gse240704/s_axis_drivers/s_axis_top_probe_drivers_annotated.csv \
  --outdir "$REANN_DIR"

echo
echo "[rescue] inspect:"
echo "  $AUDIT_DIR/bundle_download_records_audit.csv"
echo "  $INSPECT_DIR/bundle_inspection_summary.json"
echo "  $LINK_DIR/gpl23976_candidate_links.csv"
echo "  $LINK_DIR/gpl23976_candidate_payload_download_summary.json"
echo "  $PAYLOAD_INSPECT_DIR/bundle_inspection_summary.json"
echo "  $NORM_DIR/gpl23976_annotation_normalized_summary.json"
echo "  $REANN_DIR"
