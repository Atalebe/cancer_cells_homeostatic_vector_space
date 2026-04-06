#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

PYTHONPATH=. python src/cancer/03_build_gse161895_master_table.py
PYTHONPATH=. python src/cancer/03b_inspect_gse161895_realization.py
PYTHONPATH=. python src/cancer/04_define_gse161895_hrsm_proxies.py
PYTHONPATH=. python src/cancer/05_build_gse161895_state_table.py
PYTHONPATH=. python src/cancer/06_plot_gse161895_hrsm_geometry.py
PYTHONPATH=. python src/cancer/08_gse161895_unsupervised_state_domains.py

echo
echo "[done] GSE161895 branch realized through first-pass geometry/domain stage"
find results/gse161895 -maxdepth 3 -type f | sort
