#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml
from scipy.stats import fisher_exact


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


GENE_SETS = {
    "apoptosis_seed": [
        "BAX", "BAK1", "BBC3", "PMAIP1", "BCL2", "BCL2L11", "CASP3", "CASP8", "CASP9",
        "FAS", "FADD", "TNFRSF10B", "DIABLO", "APAF1", "JUN", "DDIT3"
    ],
    "hypoxia_seed": [
        "HIF1A", "EPAS1", "VEGFA", "CA9", "LDHA", "SLC2A1", "BNIP3", "PGK1", "ENO1",
        "ALDOA", "PDK1", "NDRG1"
    ],
    "erythrocyte_development_seed": [
        "HBB", "HBA1", "HBA2", "HBM", "HBQ1", "ALAS2", "AHSP", "GYPA", "GYPB",
        "KLF1", "EPOR", "SLC4A1", "TFRC", "HEMGN"
    ],
}


def fisher_from_membership(universe: set[str], target: set[str], gene_set: set[str]) -> tuple[list[list[int]], float, float]:
    a = len(target & gene_set)
    b = len(target - gene_set)
    c = len((universe - target) & gene_set)
    d = len((universe - target) - gene_set)

    table = [[a, b], [c, d]]
    odds_ratio, p_value = fisher_exact(table, alternative="greater")
    return table, odds_ratio, p_value


def main() -> None:
    cfg = read_config()
    in_dir = REPO_ROOT / cfg["results_dir"] / "d1_d2_marker_annotation"
    clean_universe_dir = REPO_ROOT / cfg["results_dir"] / "d1_d2_marker_scan_clean"
    out_dir = REPO_ROOT / cfg["results_dir"] / "d1_high_candidate_enrichment"
    out_dir.mkdir(parents=True, exist_ok=True)

    d1_high = pd.read_csv(in_dir / "d1_high_top200_clean_annotated.csv")
    clean_universe = pd.read_csv(clean_universe_dir / "d1_vs_d2_marker_scan_full_clean.csv")

    # annotated symbols
    target_symbols = set(
        d1_high["gene_symbol"].fillna("").astype(str).str.upper().replace("", pd.NA).dropna().tolist()
    )

    # need universe symbols too, via annotation table already merged into the annotated high list only
    # so derive from available gene_symbol in d1_high annotated + raw clean universe gene_ids is insufficient
    # use annotation file directly
    ann = pd.read_csv(REPO_ROOT / "data" / "reference" / "annotations" / "ensembl_gene_annotation.tsv", sep="\t")
    ann["gene_id_stripped"] = ann["gene_id"].astype(str).str.split(".", n=1).str[0]
    clean_universe["gene_id_stripped"] = clean_universe["gene_id"].astype(str).str.split(".", n=1).str[0]
    clean_universe = clean_universe.merge(
        ann[["gene_id_stripped", "gene_symbol"]],
        on="gene_id_stripped",
        how="left",
    )

    universe_symbols = set(
        clean_universe["gene_symbol"].fillna("").astype(str).str.upper().replace("", pd.NA).dropna().tolist()
    )

    rows = []
    hit_rows = []

    for set_name, genes in GENE_SETS.items():
        gene_set = {g.upper() for g in genes}
        table, odds_ratio, p_value = fisher_from_membership(universe_symbols, target_symbols, gene_set)
        hits = sorted(target_symbols & gene_set)

        rows.append({
            "gene_set": set_name,
            "n_hits_in_d1_high_top200": int(len(hits)),
            "hits": "; ".join(hits),
            "odds_ratio": float(odds_ratio) if odds_ratio is not None else None,
            "p_value": float(p_value),
            "target_size": int(len(target_symbols)),
            "gene_set_size_in_universe": int(len(universe_symbols & gene_set)),
            "contingency_a_target_and_set": int(table[0][0]),
            "contingency_b_target_not_set": int(table[0][1]),
            "contingency_c_background_set": int(table[1][0]),
            "contingency_d_background_not_set": int(table[1][1]),
        })

        for hit in hits:
            hit_rows.append({"gene_set": set_name, "gene_symbol": hit})

    summary_df = pd.DataFrame(rows).sort_values(["p_value", "n_hits_in_d1_high_top200"], ascending=[True, False])
    hits_df = pd.DataFrame(hit_rows)

    summary_df.to_csv(out_dir / "d1_high_candidate_enrichment_summary.csv", index=False)
    hits_df.to_csv(out_dir / "d1_high_candidate_enrichment_hits.csv", index=False)

    with open(out_dir / "d1_high_candidate_enrichment_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary_df.to_dict(orient="records"), f, indent=2)

    print(summary_df)


if __name__ == "__main__":
    main()
