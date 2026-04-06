#!/usr/bin/env bash
set -euo pipefail

CONFIG="configs/cancer/gse240704.yaml"

echo "[stage2] plotting HRSM geometry"
PYTHONPATH=src/cancer:. python src/cancer/06_plot_gse240704_hrsm_geometry.py "$CONFIG"

echo "[stage2] building sample annotations"
PYTHONPATH=src/cancer:. python src/cancer/07_build_gse240704_sample_annotations.py "$CONFIG"

echo "[stage2] building unsupervised state domains"
PYTHONPATH=src/cancer:. python src/cancer/08_gse240704_unsupervised_state_domains.py "$CONFIG"

echo "[stage2] building candidate irreversible shell table"
PYTHONPATH=src/cancer:. python src/cancer/09_gse240704_candidate_irreversible_core.py "$CONFIG"

echo "[stage2] running condition axis tests"
PYTHONPATH=src/cancer:. python src/cancer/10_gse240704_condition_axis_tests.py "$CONFIG"

echo "[done] stage2 finished"
