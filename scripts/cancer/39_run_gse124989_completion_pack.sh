#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

echo "============================================================"
echo "[1/5] branch threshold robustness"
PYTHONPATH=. python scripts/cancer/34_branch_threshold_robustness.py \
  --config configs/cancer/gse124989_stem_branch.yaml

echo "============================================================"
echo "[2/5] refined S k-neighbor robustness"
PYTHONPATH=. python scripts/cancer/35_s_refined_k_robustness.py \
  --config configs/cancer/gse124989_stem_branch.yaml \
  --k-values 8 10 12 15

echo "============================================================"
echo "[3/5] candidate core score size robustness"
PYTHONPATH=. python scripts/cancer/36_core_score_size_robustness.py \
  --config configs/cancer/gse124989_stem_branch.yaml \
  --min-lists 3 \
  --sizes 12 16 20

echo "============================================================"
echo "[4/5] leave-one-gene-out candidate core robustness"
PYTHONPATH=. python scripts/cancer/37_leave_one_gene_out_core_score.py \
  --config configs/cancer/gse124989_stem_branch.yaml \
  --min-lists 3 \
  --top-n 16

echo "============================================================"
echo "[5/5] sector enrichment permutation null"
PYTHONPATH=. python scripts/cancer/38_sector_enrichment_permutation_null.py \
  --config configs/cancer/gse124989_stem_branch.yaml \
  --n-perm 200 \
  --min-lists 3 \
  --top-n 16 \
  --sector malignant_reservoir

echo "============================================================"
echo "[done] gse124989 completion pack finished"
