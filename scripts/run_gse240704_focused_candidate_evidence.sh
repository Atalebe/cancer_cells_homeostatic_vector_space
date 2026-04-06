#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

PYTHONPATH=. python src/cancer/81_focused_candidate_evidence_and_regulatory_overlap.py

echo
echo "[done] focused candidate evidence follow up finished"
find results/gse240704/focused_candidate_evidence -maxdepth 3 -type f | sort
