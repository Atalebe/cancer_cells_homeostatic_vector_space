#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

echo "[reannot] rebuilding manifest-derived normalized annotation table"
PYTHONPATH=. python src/cancer/59_build_gpl23976_annotation_from_illumina_manifest.py \
  --input-file data/raw/annotations/manual_illumina/infinium-methylationepic-v-1-0-b5-manifest-file.csv \
  --outcsv results/gse240704/manual_manifest_annotation/gpl23976_annotation_from_manifest.csv \
  --outparquet results/gse240704/manual_manifest_annotation/gpl23976_annotation_from_manifest.parquet \
  --summary-json results/gse240704/manual_manifest_annotation/gpl23976_annotation_from_manifest_summary.json \
  --head-csv results/gse240704/manual_manifest_annotation/gpl23976_annotation_from_manifest_head50.csv

echo "[reannot] rerunning probe annotations with manifest-derived reference"
PYTHONPATH=. python src/cancer/53_rerun_probe_annotations_with_normalized_reference.py \
  --annotation-parquet results/gse240704/manual_manifest_annotation/gpl23976_annotation_from_manifest.parquet \
  --contrast-csvs \
    results/gse240704/d2_d3_probe_signatures/D3_vs_not_D3_top_200_probes.csv \
    results/gse240704/d2_d3_probe_signatures/D3_vs_D2_top_200_probes.csv \
    results/gse240704/d2_d3_probe_signatures/D2_vs_D1_top_200_probes.csv \
  --s-axis-driver-csv results/gse240704/s_axis_drivers/s_axis_top_probe_drivers_annotated.csv \
  --outdir results/gse240704/manifest_reannotated

echo
echo "[reannot] inspect:"
echo "  results/gse240704/manual_manifest_annotation/gpl23976_annotation_from_manifest_summary.json"
echo "  results/gse240704/manifest_reannotated/reannotation_batch_summary.json"
echo "  results/gse240704/manifest_reannotated/D3_vs_not_D3_top_200_probes_normalized_annotation.csv"
echo "  results/gse240704/manifest_reannotated/s_axis_top_probe_drivers_annotated_normalized_annotation.csv"
