#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

configs=(
  "configs/cancer/gse161895_stem_branch.yaml"
  "configs/cancer/gse240704_stem_branch.yaml"
  "configs/cancer/gse184880_recoverability_branch.yaml"
  "configs/cancer/gse161299_resistance_branch.yaml"
)

for cfg in "${configs[@]}"; do
  echo "============================================================"
  echo "[info] downloading with config: $cfg"
  PYTHONPATH=. python scripts/downloads/10_download_dataset.py \
    --config "$cfg"
done
