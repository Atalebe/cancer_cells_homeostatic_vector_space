# scripts/run_gpl23976_manifest_preinspect.sh
#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

OUTDIR="results/gse240704/manual_manifest_preinspect"
mkdir -p "${OUTDIR}"

INPUTS=(
  "data/raw/annotations/manual_illumina/infinium-methylationepic-v-1-0-b5-manifest-file.csv"
  "data/raw/annotations/manual_illumina/Infinium_MethylationEPIC_v1_0_B5_Manifest.csv.zip.zip"
  "data/raw/annotations/manual_illumina/infinium-methylationepic-v-1-0-b5-manifest-file-bpm.zip"
  "data/raw/annotations/manual_illumina/infinium-methylationepic-v1-0-missing-legacy-cpg-b3-vs-b2-annotations.zip"
)

PYTHONPATH=. python src/cancer/60_inspect_manual_illumina_manifest_inputs.py \
  --inputs "${INPUTS[@]}" \
  --outdir "${OUTDIR}"

echo
echo "[peek plain csv header area]"
PYTHONPATH=. python src/cancer/61_peek_manifest_header_region.py \
  --input "data/raw/annotations/manual_illumina/infinium-methylationepic-v-1-0-b5-manifest-file.csv" \
  --start 0 \
  --n 120

echo
echo "[inspect outputs]"
echo "  ${OUTDIR}/manual_manifest_input_inventory.csv"
echo "  ${OUTDIR}/manual_manifest_input_inspection.json"
