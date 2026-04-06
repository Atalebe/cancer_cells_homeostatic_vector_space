#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

PYTHONPATH=. python src/cancer/20_gse161895_d2_marker_scan.py
PYTHONPATH=. python src/cancer/21_gse161895_annotate_d2_marker_tables.py
PYTHONPATH=. python src/cancer/22_gse161895_d2_program_overlay_tests.py

echo
echo "[done] GSE161895 D2 marker followup finished"
find results/gse161895/d2_marker_scan -maxdepth 2 -type f | sort
find results/gse161895/d2_marker_annotation -maxdepth 2 -type f | sort
find results/gse161895/d2_program_overlay_tests -maxdepth 2 -type f | sort
