#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

PYTHONPATH=. python src/cancer/77_sim2_compactness_and_targets.py

echo
echo "[done] SIM2 follow up finished"
find results/gse240704/sim2_followup -maxdepth 3 -type f | sort
