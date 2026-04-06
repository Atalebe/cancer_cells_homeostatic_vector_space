#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

PYTHONPATH=. python src/cancer/34_gse161895_transition_trajectory_pass.py
PYTHONPATH=. python src/cancer/35_gse161895_d2_1_pseudobulk_treatment_robustness.py
PYTHONPATH=. python src/cancer/36_gse161895_reversion_candidate_scan.py

echo
echo "[done] GSE161895 transition, pseudobulk, and reversion followup finished"
find results/gse161895/transition_trajectory_pass -maxdepth 1 -type f | sort || true
find results/gse161895/d2_1_pseudobulk_treatment_robustness -maxdepth 1 -type f | sort || true
find results/gse161895/reversion_candidate_scan -maxdepth 1 -type f | sort || true
