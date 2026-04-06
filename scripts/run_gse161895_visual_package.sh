#!/usr/bin/env bash
set -euo pipefail
cd "${1:-$(pwd)}"
PYTHONPATH=. python gse161895_visual_package/src/plot_master_reversion_index.py
PYTHONPATH=. python gse161895_visual_package/src/plot_trajectory_visuals.py
PYTHONPATH=. python gse161895_visual_package/src/plot_program_shift_panels.py
PYTHONPATH=. python gse161895_visual_package/src/plot_go_mechanism_panels.py
PYTHONPATH=. python gse161895_visual_package/src/plot_variance_age_plasticity.py
echo
echo "[done] GSE161895 visual package generated under results/gse161895/visual_package/"
