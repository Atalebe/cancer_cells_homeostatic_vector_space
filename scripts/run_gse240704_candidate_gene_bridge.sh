#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

PYTHONPATH=. python src/cancer/78_write_sim2_seed_table_tex.py
PYTHONPATH=. python src/cancer/79_directional_candidate_gene_check.py

echo
echo "[done] candidate-gene bridge follow up finished"
find results/gse240704/candidate_gene_directional_check -maxdepth 3 -type f | sort
find results/gse240704/sim2_followup/tables_tex -maxdepth 2 -type f | sort
