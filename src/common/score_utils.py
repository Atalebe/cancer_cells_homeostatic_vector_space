from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd


def robust_zscore(values: pd.Series) -> pd.Series:
    x = pd.Series(values, copy=True).astype(float)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    if mad == 0 or np.isnan(mad):
        return pd.Series(np.zeros(len(x)), index=x.index)
    return (x - med) / (1.4826 * mad)


def strip_ensembl_version(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(r"\.\d+$", "", regex=True)


def build_gene_symbol_index(
    gene_map: pd.DataFrame,
    ensembl_col: str = "ensembl_id",
    symbol_col: str = "gene_symbol",
) -> pd.DataFrame:
    out = gene_map[[ensembl_col, symbol_col]].copy()
    out[ensembl_col] = strip_ensembl_version(out[ensembl_col])
    out[symbol_col] = out[symbol_col].astype(str)
    out = out.dropna().drop_duplicates()
    return out


def expression_wide_to_cell_by_gene(
    expr_wide: pd.DataFrame,
    gene_col: str = "ensembl_id",
) -> pd.DataFrame:
    cell_cols = [c for c in expr_wide.columns if c != gene_col]
    x = expr_wide.set_index(gene_col)[cell_cols].T
    x.index.name = "cell_id"
    return x


def attach_gene_symbols_to_cell_by_gene(
    x_cell_gene: pd.DataFrame,
    gene_map: pd.DataFrame,
    ensembl_col: str = "ensembl_id",
    symbol_col: str = "gene_symbol",
) -> pd.DataFrame:
    gene_index = pd.DataFrame({ensembl_col: x_cell_gene.columns})
    gene_index[ensembl_col] = strip_ensembl_version(gene_index[ensembl_col])

    gm = build_gene_symbol_index(gene_map, ensembl_col=ensembl_col, symbol_col=symbol_col)
    merged = gene_index.merge(gm, on=ensembl_col, how="left")

    new_cols = []
    seen = {}
    for ens, sym in zip(gene_index[ensembl_col], merged[symbol_col]):
        if pd.isna(sym) or str(sym).strip() == "":
            label = ens
        else:
            label = str(sym).upper().strip()
        seen[label] = seen.get(label, 0) + 1
        if seen[label] > 1:
            label = f"{label}__dup{seen[label]}"
        new_cols.append(label)

    out = x_cell_gene.copy()
    out.columns = new_cols
    return out


def score_gene_set_mean_z(
    x_cell_gene: pd.DataFrame,
    genes: Iterable[str],
) -> pd.Series:
    wanted = {str(g).upper().strip() for g in genes}
    matched = [c for c in x_cell_gene.columns if c.split("__dup")[0] in wanted]
    if len(matched) == 0:
        return pd.Series(np.zeros(x_cell_gene.shape[0]), index=x_cell_gene.index)

    sub = x_cell_gene[matched].copy()
    sub = np.log1p(sub)
    sub = sub.apply(robust_zscore, axis=0)
    return sub.mean(axis=1)


def gene_set_match_report(
    x_cell_gene: pd.DataFrame,
    genes: Iterable[str],
) -> dict:
    wanted = [str(g).upper().strip() for g in genes]
    wanted_set = set(wanted)
    available = {c.split("__dup")[0] for c in x_cell_gene.columns}
    matched = sorted(wanted_set & available)
    missing = sorted(wanted_set - available)
    return {
        "n_requested": len(wanted_set),
        "n_matched": len(matched),
        "n_missing": len(missing),
        "matched_genes": matched,
        "missing_genes": missing,
    }
