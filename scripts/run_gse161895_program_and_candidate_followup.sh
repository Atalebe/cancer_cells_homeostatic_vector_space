#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

PYTHONPATH=. python src/cancer/16_gse161895_program_overlay_tests.py
PYTHONPATH=. python src/cancer/17_gse161895_d1_high_candidate_enrichment.py

echo
echo "[done] GSE161895 program overlay and D1-high candidate enrichment finished"
find results/gse161895/program_overlay_tests -maxdepth 2 -type f | sort
find results/gse161895/d1_high_candidate_enrichment -maxdepth 2 -type f | sort
