#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

PYTHONPATH=. python src/cancer/18_gse161895_qc_artifact_challenge.py
PYTHONPATH=. python src/cancer/19_gse161895_d2_only_reanalysis.py

echo
echo "[done] GSE161895 QC artifact challenge and D2-only reanalysis finished"
find results/gse161895/qc_artifact_challenge -maxdepth 2 -type f | sort
find results/gse161895/d2_only_reanalysis -maxdepth 2 -type f | sort
