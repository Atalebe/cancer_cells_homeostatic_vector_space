#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

PYTHONPATH=. python src/cancer/25_gse161895_patient_treatment_coverage_audit.py
PYTHONPATH=. python src/cancer/26_gse161895_patient3_d2_1_treatment_sensitivity.py
PYTHONPATH=. python src/cancer/27_gse161895_hrsm_profile_visual_compare.py

echo
echo "[done] patient-treatment audit, Patient3 sensitivity, and HRSM visual comparison finished"
find results/gse161895/patient_treatment_coverage_audit -maxdepth 2 -type f | sort
find results/gse161895/patient3_d2_1_treatment_sensitivity -maxdepth 2 -type f | sort
find results/gse161895/hrsm_profile_visual_compare -maxdepth 3 -type f | sort
