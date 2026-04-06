#!/usr/bin/env bash
set -euo pipefail

CONFIG="configs/cancer/gse240704.yaml"

echo "[stage1] building metadata registry"
PYTHONPATH=src/cancer:. python src/cancer/01_build_gse240704_metadata_registry.py "$CONFIG"

echo "[stage1] auditing matrix"
PYTHONPATH=src/cancer:. python src/cancer/02_audit_gse240704_methylation_matrix.py "$CONFIG"

echo "[stage1] building beta matrix from normalized input"
PYTHONPATH=src/cancer:. python src/cancer/03_build_gse240704_beta_matrix.py "$CONFIG"

echo "[stage1] defining HRSM proxies"
PYTHONPATH=src/cancer:. python src/cancer/04_define_gse240704_hrsm_proxies.py "$CONFIG"

echo "[stage1] building state table"
PYTHONPATH=src/cancer:. python src/cancer/05_build_gse240704_state_table.py "$CONFIG"

echo "[done] stage1 finished"
