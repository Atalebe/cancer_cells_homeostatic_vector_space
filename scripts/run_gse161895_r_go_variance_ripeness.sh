#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

PYTHONPATH=. python src/cancer/28_gse161895_d2_1_r_treatment_bridge_scan.py
PYTHONPATH=. python src/cancer/29_gse161895_d2_1_r_bridge_go_enrichment.py
PYTHONPATH=. python src/cancer/30_gse161895_d2_1_variance_scaling.py
PYTHONPATH=. python src/cancer/31_gse161895_d2_1_ripeness_proxy.py

echo
echo "[done] GSE161895 R-bridge, GO, variance scaling, and ripeness followup finished"
find results/gse161895/d2_1_r_treatment_bridge -maxdepth 2 -type f | sort
find results/gse161895/d2_1_r_go_enrichment -maxdepth 2 -type f | sort
find results/gse161895/d2_1_variance_scaling -maxdepth 2 -type f | sort
find results/gse161895/d2_1_ripeness_proxy -maxdepth 2 -type f | sort
