#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

PYTHONPATH=. python src/cancer/23_gse161895_d2_immune_identity_overlay.py
PYTHONPATH=. python src/cancer/24_gse161895_d2_1_treatment_response_analysis.py

echo
echo "[done] GSE161895 D2 immune identity and D2_1 treatment followup finished"
find results/gse161895/d2_immune_identity_overlay -maxdepth 2 -type f | sort
find results/gse161895/d2_1_treatment_response -maxdepth 2 -type f | sort
