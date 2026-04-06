#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

PYTHONPATH=. python src/cancer/75_d3_21q_locus_check.py
PYTHONPATH=. python src/cancer/76_plot_d3_21q_loci.py || true

echo
echo "[done] 21q locus follow up finished"
echo "[check] results/gse240704/chr21_locus_followup/"
find results/gse240704/chr21_locus_followup -maxdepth 2 -type f | sort
