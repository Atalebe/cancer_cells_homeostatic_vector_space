#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

PYTHONPATH=. python src/cancer/32_gse161895_plasticity_reversion_diagnostics.py
PYTHONPATH=. python src/cancer/33_gse161895_age_axis_proxy.py

echo
echo "[done] GSE161895 plasticity and age-axis followup finished"
echo "results/gse161895/plasticity_reversion_diagnostics/plasticity_reversion_summary.json"
echo "results/gse161895/plasticity_reversion_diagnostics/plasticity_group_summary.csv"
echo "results/gse161895/plasticity_reversion_diagnostics/plasticity_reversion_tests.csv"
echo "results/gse161895/age_axis_proxy/age_axis_summary.json"
echo "results/gse161895/age_axis_proxy/age_axis_group_summary.csv"
echo "results/gse161895/age_axis_proxy/age_axis_tests.csv"
