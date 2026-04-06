from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import wilcoxon


REPO = Path(__file__).resolve().parents[2]
PROC = REPO / "data" / "processed" / "gse161895"
D2_DIR = REPO / "results" / "gse161895" / "d2_only_reanalysis"
OUTDIR = REPO / "results" / "gse161895" / "d2_1_pseudobulk_treatment_robustness"


def load_counts() -> tuple[sparse.csr_matrix, list[str], list[str]]:
    npz_path = PROC / "counts_sparse.npz"
    gene_path = PROC / "gene_ids.csv"
    cell_path = PROC / "cell_ids.csv"

    if not npz_path.exists():
        raise FileNotFoundError(f"Missing sparse count matrix: {npz_path}")
    if not gene_path.exists():
        raise FileNotFoundError(f"Missing gene id file: {gene_path}")
    if not cell_path.exists():
        raise FileNotFoundError(f"Missing cell id file: {cell_path}")

    counts = sparse.load_npz(npz_path).tocsr()

    gene_df = pd.read_csv(gene_path)
    cell_df = pd.read_csv(cell_path)

    gene_col = gene_df.columns[0]
    cell_col = cell_df.columns[0]

    genes = gene_df[gene_col].astype(str).tolist()
    cells = cell_df[cell_col].astype(str).tolist()

    if counts.shape[0] != len(cells):
        if counts.shape[1] == len(cells):
            counts = counts.T.tocsr()
        else:
            raise ValueError(
                f"Could not align sparse matrix with cells. "
                f"counts shape={counts.shape}, n_cells={len(cells)}"
            )

    if counts.shape[1] != len(genes):
        raise ValueError(
            f"Could not align sparse matrix with genes. "
            f"counts shape={counts.shape}, n_genes={len(genes)}"
        )

    return counts, cells, genes


def load_metadata() -> pd.DataFrame:
    meta = pd.read_parquet(PROC / "cell_metadata_registry_aligned.parquet")
    d2 = pd.read_parquet(D2_DIR / "d2_state_domains.parquet")

    cols = [c for c in ["cell_id", "patient_or_donor_id", "source_class", "treatment_state"] if c in meta.columns]
    meta = meta[cols].copy()
    meta["cell_id"] = meta["cell_id"].astype(str)

    d2["cell_id"] = d2["cell_id"].astype(str)
    d2 = d2[["cell_id", "d2_subdomain"]].copy()

    merged = meta.merge(d2, on="cell_id", how="left")
    merged["treatment_state"] = merged["treatment_state"].fillna("").astype(str).str.lower().str.strip()
    merged["patient_or_donor_id"] = merged["patient_or_donor_id"].fillna("").astype(str)

    return merged


def compute_library_normalized_log1p(sub_counts: sparse.csr_matrix) -> np.ndarray:
    lib = np.asarray(sub_counts.sum(axis=1)).ravel().astype(float)
    lib[lib <= 0] = 1.0
    scale = 1e4 / lib
    norm = sub_counts.multiply(scale[:, None])
    return np.log1p(norm.toarray())


def build_program_sets() -> dict[str, list[str]]:
    return {
        "antigen_presentation": ["B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5"],
        "immune_leaning": ["IGHM", "MS4A1", "FCRL3", "BTLA", "KLRB1", "FOXO1", "FCMR", "SKAP1", "CRTAM", "THEMIS", "GZMM", "P2RY10", "MAL"],
        "malignant_output_core": ["ACTB", "ACTG1", "EEF1A1", "DDX5", "HSPA8", "RACK1", "HSP90AB1", "ENO1", "FTL", "GAPDH", "RPL3", "RPL4", "RPLP0", "RPS6"],
        "mitochondrial_encoded": ["MT-RNR1", "MT-RNR2", "MT-CO1", "MT-CO2", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND4"],
        "ribosomal_translation": ["RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL13", "RPL29", "RPS3", "RPS4X", "RPS6", "RPS7", "RPS23"],
    }


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    counts, cell_ids, gene_ids = load_counts()
    meta = load_metadata()

    cell_index = pd.Index(cell_ids, name="cell_id")
    gene_index = pd.Index(gene_ids, name="gene_id")

    meta = meta[meta["cell_id"].isin(cell_index)].copy()
    meta = meta[meta["d2_subdomain"] == "D2_1"].copy()
    meta = meta[meta["treatment_state"].isin(["treated", "untreated"])].copy()

    patient_treatment = pd.crosstab(meta["patient_or_donor_id"], meta["treatment_state"])
    paired_patients = patient_treatment.index[
        (patient_treatment.get("treated", 0) > 0) & (patient_treatment.get("untreated", 0) > 0)
    ].tolist()

    if not paired_patients:
        raise ValueError("No patients with both treated and untreated D2_1 cells were found.")

    meta = meta[meta["patient_or_donor_id"].isin(paired_patients)].copy()

    cell_pos = cell_index.get_indexer(meta["cell_id"])
    if np.any(cell_pos < 0):
        raise ValueError("Some D2_1 cells could not be found in the sparse matrix.")

    sub_counts = counts[cell_pos, :]
    expr = compute_library_normalized_log1p(sub_counts)

    expr_df = pd.DataFrame(expr, index=meta["cell_id"].tolist(), columns=gene_index)

    patient_gene_rows = []
    for patient in paired_patients:
        for treatment in ["treated", "untreated"]:
            ids = meta.loc[
                (meta["patient_or_donor_id"] == patient) &
                (meta["treatment_state"] == treatment),
                "cell_id"
            ].tolist()
            if not ids:
                continue
            mean_expr = expr_df.loc[ids].mean(axis=0)
            row = pd.DataFrame({
                "patient_or_donor_id": patient,
                "treatment_state": treatment,
                "gene_id": gene_index,
                "pseudobulk_mean_log1p": mean_expr.values,
            })
            patient_gene_rows.append(row)

    gene_long = pd.concat(patient_gene_rows, ignore_index=True)

    treated = gene_long[gene_long["treatment_state"] == "treated"].rename(
        columns={"pseudobulk_mean_log1p": "treated_mean_log1p"}
    )
    untreated = gene_long[gene_long["treatment_state"] == "untreated"].rename(
        columns={"pseudobulk_mean_log1p": "untreated_mean_log1p"}
    )

    paired_gene = treated.merge(
        untreated,
        on=["patient_or_donor_id", "gene_id"],
        how="inner",
    )
    paired_gene["delta_treated_minus_untreated"] = (
        paired_gene["treated_mean_log1p"] - paired_gene["untreated_mean_log1p"]
    )

    gene_summary = (
        paired_gene.groupby("gene_id")
        .agg(
            n_patients=("patient_or_donor_id", "nunique"),
            mean_delta=("delta_treated_minus_untreated", "mean"),
            median_delta=("delta_treated_minus_untreated", "median"),
        )
        .reset_index()
    )

    pvals = []
    for gene_id, sub in paired_gene.groupby("gene_id", sort=False):
        vals = sub["delta_treated_minus_untreated"].to_numpy(dtype=float)
        if len(vals) < 2 or np.allclose(vals, 0):
            pvals.append((gene_id, np.nan))
            continue
        try:
            p = wilcoxon(vals, zero_method="wilcox", alternative="two-sided").pvalue
        except ValueError:
            p = np.nan
        pvals.append((gene_id, p))

    pval_df = pd.DataFrame(pvals, columns=["gene_id", "p_value"])
    gene_summary = gene_summary.merge(pval_df, on="gene_id", how="left")
    gene_summary["abs_median_delta"] = gene_summary["median_delta"].abs()
    gene_summary = gene_summary.sort_values(["p_value", "abs_median_delta"], ascending=[True, False])

    gene_summary.to_csv(OUTDIR / "d2_1_pseudobulk_gene_robustness.csv", index=False)

    programs = build_program_sets()
    gene_symbol_map = {}
    ann_path = REPO / "data" / "reference" / "annotations" / "ensembl_gene_annotation.tsv"
    if ann_path.exists():
        ann = pd.read_csv(ann_path, sep="\t")
        if "gene_id" in ann.columns and "gene_symbol" in ann.columns:
            ann["gene_id"] = ann["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
            gene_symbol_map = ann.drop_duplicates("gene_id").set_index("gene_id")["gene_symbol"].to_dict()

    stripped_gene_ids = pd.Series(gene_index.astype(str)).str.replace(r"\.\d+$", "", regex=True).tolist()
    stripped_to_full = dict(zip(stripped_gene_ids, gene_index.tolist()))

    program_rows = []
    for program_name, symbols in programs.items():
        program_gene_ids = []
        for gid in gene_summary["gene_id"].astype(str):
            stripped = gid.replace(".0", "")
            sym = gene_symbol_map.get(stripped, "")
            if sym in symbols:
                program_gene_ids.append(gid)
        if not program_gene_ids:
            continue

        sub = paired_gene[paired_gene["gene_id"].isin(program_gene_ids)].copy()
        per_patient = (
            sub.groupby("patient_or_donor_id")
            .agg(program_delta_mean=("delta_treated_minus_untreated", "mean"))
            .reset_index()
        )
        vals = per_patient["program_delta_mean"].to_numpy(dtype=float)
        if len(vals) >= 2 and not np.allclose(vals, 0):
            try:
                p = wilcoxon(vals, zero_method="wilcox", alternative="two-sided").pvalue
            except ValueError:
                p = np.nan
        else:
            p = np.nan

        program_rows.append({
            "program": program_name,
            "n_genes": len(program_gene_ids),
            "n_patients": len(per_patient),
            "median_patient_delta": float(np.median(vals)) if len(vals) else np.nan,
            "mean_patient_delta": float(np.mean(vals)) if len(vals) else np.nan,
            "p_value": p,
            "genes": "; ".join(sorted({
                gene_symbol_map.get(g.replace(".0", "").split(".")[0], g) for g in program_gene_ids
            })),
        })

    program_df = pd.DataFrame(program_rows).sort_values("p_value", ascending=True, na_position="last")
    program_df.to_csv(OUTDIR / "d2_1_pseudobulk_program_robustness.csv", index=False)

    summary = {
        "paired_patients": paired_patients,
        "n_paired_patients": len(paired_patients),
        "n_d2_1_cells_used": int(meta.shape[0]),
        "outputs": {
            "gene_csv": str(OUTDIR / "d2_1_pseudobulk_gene_robustness.csv"),
            "program_csv": str(OUTDIR / "d2_1_pseudobulk_program_robustness.csv"),
        },
    }
    with open(OUTDIR / "d2_1_pseudobulk_robustness_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
