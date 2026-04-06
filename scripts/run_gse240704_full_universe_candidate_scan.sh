#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

PYTHONPATH=. python src/cancer/80_full_universe_candidate_scan.py

echo
echo "[done] full-universe candidate scan finished"
find results/gse240704/full_universe_candidate_scan -maxdepth 3 -type f | sort
