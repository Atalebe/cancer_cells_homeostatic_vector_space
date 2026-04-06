#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

PYTHONPATH=. python src/cancer/83_ptch1_olig2_tfbs_closing_analysis.py

echo
echo "[done] PTCH1 and OLIG2 TFBS closing analysis finished"
find results/gse240704/ptch1_olig2_tfbs_closing_analysis -maxdepth 3 -type f | sort
