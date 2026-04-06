#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

PYTHONPATH=. python src/cancer/82_matched_background_regulatory_enrichment.py

echo
echo "[done] matched-background regulatory enrichment finished"
find results/gse240704/matched_background_regulatory_enrichment -maxdepth 3 -type f | sort
