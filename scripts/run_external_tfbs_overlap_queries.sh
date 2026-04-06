#!/usr/bin/env bash
set -euo pipefail

cd ~/cancer_cells_hrsm_sims

INDIR="data/reference/external_interval_queries"
OUTDIR="results/gse240704/external_tfbs_overlap"
mkdir -p "$OUTDIR"/{remap_json,ucsc_json,logs}

INTERVALS="${INDIR}/gse240704_closing_intervals_hg19.tsv"

# ReMap version/datatype may need adjusting if the service changes names.
# The REST docs confirm the endpoint shape, but not the enum values in the snippet.
# This script starts with a widely plausible default and records failures cleanly.
REMAP_VERSION="${REMAP_VERSION:-2022}"
REMAP_ASSEMBLY="${REMAP_ASSEMBLY:-hg19}"
REMAP_DATATYPE="${REMAP_DATATYPE:-ChIP-seq}"

tail -n +2 "$INTERVALS" | while IFS=$'\t' read -r label chr start end assembly note; do
    region="${chr}:${start}-${end}"

    echo "[info] querying ReMap for ${label} ${region}" | tee -a "$OUTDIR/logs/remap.log"
    remap_url="https://remap-rest.univ-amu.fr/api/V1/get_peaks/${REMAP_VERSION}/${REMAP_ASSEMBLY}/${REMAP_DATATYPE}/${region}.json"
    curl -L --fail --silent --show-error "$remap_url" \
        -o "$OUTDIR/remap_json/${label}.json" \
        || echo "[warn] ReMap query failed for ${label} ${region}" | tee -a "$OUTDIR/logs/remap.log"

    echo "[info] querying UCSC hg19 ENCODE TFBS clustered track for ${label} ${region}" | tee -a "$OUTDIR/logs/ucsc.log"
    ucsc_url="https://api.genome.ucsc.edu/getData/track?genome=hg19;track=wgEncodeRegTfbsClusteredV3;chrom=${chr};start=${start};end=${end}"
    curl -L --fail --silent --show-error "$ucsc_url" \
        -o "$OUTDIR/ucsc_json/${label}.json" \
        || echo "[warn] UCSC query failed for ${label} ${region}" | tee -a "$OUTDIR/logs/ucsc.log"
done

echo
echo "[done] external overlap queries finished"
find "$OUTDIR" -maxdepth 2 -type f | sort
