#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

MANIFEST="data/raw/annotations/manual_illumina/infinium-methylationepic-v-1-0-b5-manifest-file.csv"
OUTDIR="results/gse240704/manual_manifest_annotation"
LOGDIR="logs"

mkdir -p "$OUTDIR" "$LOGDIR"

echo "[manifest] quick file sanity"
python - <<'PY'
from pathlib import Path

p = Path("data/raw/annotations/manual_illumina/infinium-methylationepic-v-1-0-b5-manifest-file.csv")
print("exists:", p.exists())
print("size_bytes:", p.stat().st_size if p.exists() else None)

with p.open("r", encoding="utf-8", errors="replace") as f:
    for i in range(12):
        line = f.readline()
        if not line:
            break
        print(f"{i:05d}: {line[:220].rstrip()}")
PY

echo
echo "[manifest] building normalized GPL23976 table"
nice -n 10 env PYTHONPATH=. python src/cancer/59_build_gpl23976_annotation_from_illumina_manifest.py \
  --input "$MANIFEST" \
  --outcsv "$OUTDIR/gpl23976_annotation_from_manifest.csv" \
  --outparquet "$OUTDIR/gpl23976_annotation_from_manifest.parquet" \
  --summary-json "$OUTDIR/gpl23976_annotation_from_manifest_summary.json" \
  --preview-head-csv "$OUTDIR/gpl23976_annotation_from_manifest_head50.csv" \
  2>&1 | tee "$LOGDIR/gpl23976_manual_manifest_parser.log"

echo
echo "[manifest] inspect:"
echo "  $OUTDIR/gpl23976_annotation_from_manifest_summary.json"
echo "  $OUTDIR/gpl23976_annotation_from_manifest_head50.csv"
echo "  $OUTDIR/gpl23976_annotation_from_manifest.parquet"
echo "  $LOGDIR/gpl23976_manual_manifest_parser.log"
