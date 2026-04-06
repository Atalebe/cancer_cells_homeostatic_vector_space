#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

PYTHONPATH=. python src/cancer/85_external_recurrent_tf_and_mutation_check.py

echo
echo "[done] external recurrent TF and mutation check finished"
find results/gse240704/external_tf_recurrence_and_mutation_check -maxdepth 3 -type f | sort
