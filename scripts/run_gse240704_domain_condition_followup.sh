#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

OUTDIR="results/gse240704/domain_condition_followup"
PLOTDIR="${OUTDIR}/plots"

mkdir -p "$OUTDIR" "$PLOTDIR"

echo "[followup] domain-condition enrichment and shell audit"
PYTHONPATH=. python src/cancer/27_domain_vs_condition_enrichment_gse240704.py \
  --state-domains-csv results/gse240704/state_domains/state_domains.csv \
  --sample-annotations-curated data/processed/gse240704/sample_annotations_curated.parquet \
  --stable-shell-csv results/gse240704/inspection/stable_shell/stable_shell_samples_inspection.csv \
  --outdir "$OUTDIR"

echo "[followup] plotting domain and condition panels"
PYTHONPATH=. python src/cancer/28_plot_gse240704_domain_condition_panels.py \
  --state-domains-csv results/gse240704/state_domains/state_domains.csv \
  --sample-annotations-curated data/processed/gse240704/sample_annotations_curated.parquet \
  --outdir "$PLOTDIR"

echo "[followup] inspect:"
echo "  ${OUTDIR}/domain_condition_summary.json"
echo "  ${OUTDIR}/domain_mgmt_enrichment_fisher.csv"
echo "  ${OUTDIR}/domain_axis_summary.csv"
echo "  ${OUTDIR}/domain_axis_kruskal_tests.csv"
echo "  ${OUTDIR}/stable_shell_vs_rest_curated_tests.csv"
echo "  ${OUTDIR}/stable_shell_domain_counts.csv"
echo "  ${PLOTDIR}/mgmt_counts_by_state_domain.png"
