#!/usr/bin/env bash
set -euo pipefail

CONFIG="configs/cancer/gse240704.yaml"

echo "[resume] defining HRSM proxies"
PYTHONPATH=src/cancer:. python src/cancer/04_define_gse240704_hrsm_proxies.py "$CONFIG"

echo "[resume] building state table"
PYTHONPATH=src/cancer:. python src/cancer/05_build_gse240704_state_table.py "$CONFIG"

echo "[done] resume finished"
