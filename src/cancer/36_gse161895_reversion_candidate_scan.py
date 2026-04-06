from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import mannwhitneyu


REPO = Path(__file__).resolve().parents[2]
PROC = REPO / "data" / "processed" / "gse161895"
D2_DIR = REPO / "results" / "gse161895" / "d2_only_reanalysis"
OUTDIR = REPO / "results" / "gse161895" / "reversion_candidate_scan"


def load_counts() -> tuple[sparse.csr_matrix, list[str], list[str]]:
    npz_path = PROC / "counts_sparse.npz"
    gene_path = PROC / "gene_ids.csv"
    cell_path = PROC / "cell_ids.csv"

    counts = sparse.load_npz(npz_path).tocsr()

    gene_df = pd.read_csv(gene_path)
    cell_df = pd.read_csv(cell_path)

    genes = gene_df.iloc[:, 0].astype(str).tolist()
    cells = cell_df.iloc[:, 0].astype(str).tolist()

    if counts.shape[0] != len(cells):
        if counts.shape[1] == len(cells):
            counts = counts.T.tocsr()
        else:
            raise ValueError(f"Count matrix shape {counts.shape} does not align with cell ids.")
    if counts.shape[1] != len(genes):
        raise ValueError(f"Count matrix shape {counts.shape} does not align with gene ids.")

    return counts, cells, genes


def load_metadata() -> pd.DataFrame:
    meta = pd.read_parquet(PROC / "cell_metadata_registry_aligned.parquet")
    d2 = pd.read_parquet(D2_DIR / "d2_state_domains.parquet")

    meta["cell_id"] = meta["cell_id"].astype(str)
    d2["cell_id"] = d2["cell_id"].astype(str)

    keep_cols = [c for c in ["cell_id", "patient_or_donor_id", "source_class", "treatment_state"] if c in meta.columns]
    meta = meta[keep_cols].copy()
    merged = meta.merge(d2[["cell_id", "d2_subdomain"]], on="cell_id", how="left")

    merged["source_class"] = merged["source_class"].fillna("").astype(str).str.lower().str.strip()
    merged["treatment_state"] = merged["treatment_state"].fillna("").astype(str).str.lower().str.strip()
    return merged


def compute_log1p_norm(sub_counts: sparse.csr_matrix) -> np.ndarray:
    lib = np.asarray(sub_counts.sum(axis=1)).ravel().astype(float)
    lib[lib <= 0] = 1.0
    norm = sub_counts.multiply((1e4 / lib)[:, None])
    return np.log1p(norm.toarray())


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    counts, cells, genes = load_counts()
    meta = load_metadata()

    cell_index = pd.Index(cells, name="cell_id")
    gene_index = pd.Index(genes, name="gene_id")

    meta = meta[meta["cell_id"].isin(cell_index)].copy()

    donor_ids = meta.loc[meta["source_class"] == "normal_donor", "cell_id"].tolist()
    treated_ids = meta.loc[
        (meta["d2_subdomain"] == "D2_1") & (meta["treatment_state"] == "treated"),
        "cell_id"
    ].tolist()
    untreated_ids = meta.loc[
        (meta["d2_subdomain"] == "D2_1") & (meta["treatment_state"] == "untreated"),
        "cell_id"
    ].tolist()

    if not donor_ids:
        raise ValueError("No donor cells found.")
    if not treated_ids or not untreated_ids:
        raise ValueError("Need both treated and untreated D2_1 cells for reversion scan.")

    all_ids = donor_ids + treated_ids + untreated_ids
    cell_pos = cell_index.get_indexer(all_ids)
    sub_counts = counts[cell_pos, :]
    expr = compute_log1p_norm(sub_counts)

    groups = (
        ["donor"] * len(donor_ids) +
        ["treated"] * len(treated_ids) +
        ["untreated"] * len(untreated_ids)
    )

    expr_df = pd.DataFrame(expr, index=all_ids, columns=gene_index)
    group_s = pd.Series(groups, index=all_ids, name="group")

    ann_path = REPO / "data" / "reference" / "annotations" / "ensembl_gene_annotation.tsv"
    ann_map = {}
    if ann_path.exists():
        ann = pd.read_csv(ann_path, sep="\t")
        if "gene_id" in ann.columns and "gene_symbol" in ann.columns:
            ann["gene_id"] = ann["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
            ann_map = ann.drop_duplicates("gene_id").set_index("gene_id")["gene_symbol"].to_dict()

    results = []
    donor_mask = group_s == "donor"
    treated_mask = group_s == "treated"
    untreated_mask = group_s == "untreated"

    for gene_id in expr_df.columns:
        x_d = expr_df.loc[donor_mask, gene_id].to_numpy(dtype=float)
        x_t = expr_df.loc[treated_mask, gene_id].to_numpy(dtype=float)
        x_u = expr_df.loc[untreated_mask, gene_id].to_numpy(dtype=float)

        med_d = float(np.median(x_d))
        med_t = float(np.median(x_t))
        med_u = float(np.median(x_u))

        dist_treated_to_donor = abs(med_t - med_d)
        dist_untreated_to_donor = abs(med_u - med_d)
        reversion_gain = dist_untreated_to_donor - dist_treated_to_donor

        try:
            p_treat = mannwhitneyu(x_t, x_u, alternative="two-sided").pvalue
        except ValueError:
            p_treat = np.nan

        direction = "toward_donor" if reversion_gain > 0 else "away_from_donor_or_neutral"

        stripped = str(gene_id).replace(".0", "")
        stripped = stripped.split(".")[0]
        gene_symbol = ann_map.get(stripped, stripped)

        results.append({
            "gene_id": gene_id,
            "gene_symbol": gene_symbol,
            "median_donor": med_d,
            "median_treated": med_t,
            "median_untreated": med_u,
            "treated_minus_untreated": med_t - med_u,
            "dist_treated_to_donor": dist_treated_to_donor,
            "dist_untreated_to_donor": dist_untreated_to_donor,
            "reversion_gain": reversion_gain,
            "direction": direction,
            "p_treatment": p_treat,
        })

    res = pd.DataFrame(results)
    res["abs_reversion_gain"] = res["reversion_gain"].abs()
    res = res.sort_values(["direction", "reversion_gain"], ascending=[True, False])

    res.to_csv(OUTDIR / "reversion_candidate_full_scan.csv", index=False)
    res[res["reversion_gain"] > 0].sort_values("reversion_gain", ascending=False).head(200).to_csv(
        OUTDIR / "reversion_candidate_top_toward_donor.csv", index=False
    )
    res[res["reversion_gain"] <= 0].sort_values("reversion_gain", ascending=True).head(200).to_csv(
        OUTDIR / "reversion_candidate_top_away_or_neutral.csv", index=False
    )

    summary = {
        "n_donor_cells": len(donor_ids),
        "n_treated_d2_1_cells": len(treated_ids),
        "n_untreated_d2_1_cells": len(untreated_ids),
        "n_genes_scanned": int(res.shape[0]),
        "n_toward_donor_genes": int((res["reversion_gain"] > 0).sum()),
        "outputs": {
            "full_csv": str(OUTDIR / "reversion_candidate_full_scan.csv"),
            "toward_donor_csv": str(OUTDIR / "reversion_candidate_top_toward_donor.csv"),
            "away_or_neutral_csv": str(OUTDIR / "reversion_candidate_top_away_or_neutral.csv"),
        },
    }
    with open(OUTDIR / "reversion_candidate_scan_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
