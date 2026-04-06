#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

bash scripts/run_external_tfbs_overlap_queries.sh
PYTHONPATH=. python src/cancer/84_summarize_external_tfbs_overlaps.py

echo
echo "[done] external TFBS follow up complete"
find results/gse240704/external_tfbs_overlap -maxdepth 3 -type f | sort
