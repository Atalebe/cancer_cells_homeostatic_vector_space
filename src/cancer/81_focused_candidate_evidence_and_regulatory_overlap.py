#!/usr/bin/env python3

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
DATASET = "gse240704"

FULL_UNIVERSE_TABLE = (
    REPO_ROOT
    / "results"
    / DATASET
    / "d2_d3_probe_signatures"
    / "D3_vs_not_D3_all_probes.csv"
)

TOP200_DIRECTIONAL_TABLE = (
    REPO_ROOT
    / "results"
    / DATASET
    / "directional_biology_followup"
    / "D3_vs_not_D3_top_200_probes_with_direction.csv"
)

MANIFEST_TABLE = (
    REPO_ROOT
    / "results"
    / DATASET
    / "manual_manifest_annotation"
    / "gpl23976_annotation_from_manifest.parquet"
)

OUTDIR = (
    REPO_ROOT
    / "results"
    / DATASET
    / "focused_candidate_evidence"
)

PLOTDIR = OUTDIR / "plots"
TABLETEXDIR = OUTDIR / "tables_tex"
LOGBOOKDIR = OUTDIR / "logbook"

ENTRY_NUMBER = 31

# Six detected bridge genes from the full-universe scan
FOCUSED_GENES = ["PAX6", "GLI2", "BMP4", "PTCH1", "SHH", "OLIG2"]

# Core local anchor
LOCAL_ANCHOR = "SIM2"


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [
        re.sub(r"[^a-z0-9]+", "_", str(c).strip().lower()).strip("_")
        for c in out.columns
    ]
    return out


def deduplicate_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df.loc[:, ~df.columns.duplicated()].copy()


def find_first_existing(df: pd.DataFrame, candidates: list[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def split_gene_tokens(value: object) -> list[str]:
    if pd.isna(value):
        return []
    s = str(value).strip()
    if not s:
        return []
    out = []
    for part in re.split(r"[;,/|]+", s):
        token = part.strip()
        if token:
            out.append(token)
    return out


def latex_escape(text: object) -> str:
    s = "" if pd.isna(text) else str(text)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    for old, new in replacements.items():
        s = s.replace(old, new)
    return s


def load_table_with_manifest(base_table: Path) -> pd.DataFrame:
    if not base_table.exists():
        raise FileNotFoundError(f"Missing table: {base_table}")
    if not MANIFEST_TABLE.exists():
        raise FileNotFoundError(f"Missing manifest table: {MANIFEST_TABLE}")

    base = normalize_columns(pd.read_csv(base_table))
    manifest = normalize_columns(pd.read_parquet(MANIFEST_TABLE))

    base = deduplicate_columns(base)
    manifest = deduplicate_columns(manifest)

    probe_base = find_first_existing(base, ["probe_id", "ilmnid", "ilmn_id", "name", "id_ref", "id"])
    probe_ann = find_first_existing(manifest, ["probe_id", "ilmnid", "ilmn_id", "name", "id_ref", "id"])

    if probe_base is None or probe_ann is None:
        raise ValueError("Could not identify probe id columns for merge.")

    base[probe_base] = base[probe_base].astype(str).str.strip()
    manifest[probe_ann] = manifest[probe_ann].astype(str).str.strip()

    merged = base.merge(
        manifest,
        left_on=probe_base,
        right_on=probe_ann,
        how="left",
        suffixes=("_dir", "_ann"),
    )
    merged = normalize_columns(merged)
    merged = deduplicate_columns(merged)
    return merged


def choose_gene_col(df: pd.DataFrame) -> Optional[str]:
    return find_first_existing(
        df,
        [
            "ucsc_refgene_name_ann",
            "ucsc_refgene_name_dir",
            "ucsc_refgene_name",
            "gene_symbol",
            "gene_symbols",
            "genes",
        ],
    )


def choose_cols(df: pd.DataFrame) -> dict[str, Optional[str]]:
    return {
        "direction": find_first_existing(df, ["direction", "probe_direction", "contrast_direction"]),
        "chr": find_first_existing(df, ["chromosome_ann", "chromosome_dir", "chromosome", "chromosome_std", "chr"]),
        "pos": find_first_existing(df, ["mapinfo_ann", "mapinfo_dir", "mapinfo", "position_std", "position"]),
        "delta_median": find_first_existing(df, ["delta_median_a_minus_b"]),
        "delta_mean": find_first_existing(df, ["delta_mean_a_minus_b"]),
        "abs_delta_median": find_first_existing(df, ["abs_delta_median"]),
        "q_value": find_first_existing(df, ["q_value_bh"]),
        "p_value": find_first_existing(df, ["p_value"]),
        "rank_biserial": find_first_existing(df, ["rank_biserial"]),
        "refgene_group": find_first_existing(
            df,
            [
                "refgene_group_ann",
                "refgene_group_dir",
                "ucsc_refgene_group_ann",
                "ucsc_refgene_group_dir",
                "refgene_group",
            ],
        ),
        "cpg_relation": find_first_existing(
            df,
            [
                "relation_to_cpg_island_ann",
                "relation_to_cpg_island_dir",
                "relation_to_ucsc_cpg_island_ann",
                "relation_to_ucsc_cpg_island_dir",
                "relation_to_cpg_island",
            ],
        ),
        "reg_feature_group": find_first_existing(
            df,
            [
                "regulatory_feature_group_ann",
                "regulatory_feature_group_dir",
                "regulatory_feature_group",
            ],
        ),
        "reg_feature_name": find_first_existing(
            df,
            [
                "regulatory_feature_name_ann",
                "regulatory_feature_name_dir",
                "regulatory_feature_name",
            ],
        ),
        "phantom4": find_first_existing(df, ["phantom4_enhancers_ann", "phantom4_enhancers_dir", "phantom4_enhancers"]),
        "dnase": find_first_existing(df, ["dnase_name"]),
        "openchromatin": find_first_existing(df, ["openchromatin_name"]),
        "tfbs": find_first_existing(df, ["tfbs_name"]),
    }


def row_matches_gene(row: pd.Series, gene_col: str, gene: str) -> bool:
    tokens = {tok.upper() for tok in split_gene_tokens(row.get(gene_col))}
    return gene.upper() in tokens


def extract_gene_probe_rows(df: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    gene_col = choose_gene_col(df)
    if gene_col is None:
        raise ValueError("Could not identify gene annotation column.")

    cols = choose_cols(df)

    rows = []
    for _, row in df.iterrows():
        gene_tokens = {tok.upper() for tok in split_gene_tokens(row.get(gene_col))}
        overlap = sorted(gene_tokens.intersection({g.upper() for g in genes}))
        if not overlap:
            continue

        for gene in overlap:
            rows.append(
                {
                    "candidate_gene": gene,
                    "probe_id": row.get("probe_id"),
                    "direction": row.get(cols["direction"]) if cols["direction"] else None,
                    "chromosome": row.get(cols["chr"]) if cols["chr"] else None,
                    "position": row.get(cols["pos"]) if cols["pos"] else None,
                    "delta_median_a_minus_b": row.get(cols["delta_median"]) if cols["delta_median"] else None,
                    "delta_mean_a_minus_b": row.get(cols["delta_mean"]) if cols["delta_mean"] else None,
                    "abs_delta_median": row.get(cols["abs_delta_median"]) if cols["abs_delta_median"] else None,
                    "q_value_bh": row.get(cols["q_value"]) if cols["q_value"] else None,
                    "p_value": row.get(cols["p_value"]) if cols["p_value"] else None,
                    "rank_biserial": row.get(cols["rank_biserial"]) if cols["rank_biserial"] else None,
                    "refgene_group": row.get(cols["refgene_group"]) if cols["refgene_group"] else None,
                    "relation_to_cpg_island": row.get(cols["cpg_relation"]) if cols["cpg_relation"] else None,
                    "regulatory_feature_group": row.get(cols["reg_feature_group"]) if cols["reg_feature_group"] else None,
                    "regulatory_feature_name": row.get(cols["reg_feature_name"]) if cols["reg_feature_name"] else None,
                    "phantom4_enhancers": row.get(cols["phantom4"]) if cols["phantom4"] else None,
                    "dnase_name": row.get(cols["dnase"]) if cols["dnase"] else None,
                    "openchromatin_name": row.get(cols["openchromatin"]) if cols["openchromatin"] else None,
                    "tfbs_name": row.get(cols["tfbs"]) if cols["tfbs"] else None,
                    "raw_gene_annotation": row.get(gene_col),
                }
            )

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    out["_sort_abs"] = pd.to_numeric(out["abs_delta_median"], errors="coerce")
    out = out.sort_values(["candidate_gene", "_sort_abs"], ascending=[True, False]).drop(columns=["_sort_abs"])
    return out.reset_index(drop=True)


def summarize_by_gene(hits: pd.DataFrame) -> pd.DataFrame:
    if hits.empty:
        return pd.DataFrame(
            columns=[
                "candidate_gene",
                "n_probes",
                "best_abs_delta_median",
                "min_q_value_bh",
                "directions",
                "chromosomes",
                "regulatory_feature_groups",
                "tfbs_any",
                "enhancer_any",
                "dnase_any",
                "openchromatin_any",
            ]
        )

    rows = []
    for gene, sub in hits.groupby("candidate_gene", sort=True):
        best_abs = pd.to_numeric(sub["abs_delta_median"], errors="coerce").max()
        min_q = pd.to_numeric(sub["q_value_bh"], errors="coerce").min()

        directions = ",".join(sorted(set(sub["direction"].dropna().astype(str).tolist())))
        chrs = ",".join(sorted(set(sub["chromosome"].dropna().astype(str).tolist())))

        reg_groups = sorted(set(sub["regulatory_feature_group"].dropna().astype(str).tolist()))
        reg_groups_str = ",".join(reg_groups)

        tfbs_any = sub["tfbs_name"].fillna("").astype(str).str.strip().ne("").any()
        enhancer_any = (
            sub["phantom4_enhancers"].fillna("").astype(str).str.strip().ne("").any()
            or sub["regulatory_feature_name"].fillna("").astype(str).str.contains("enhancer", case=False, na=False).any()
        )
        dnase_any = sub["dnase_name"].fillna("").astype(str).str.strip().ne("").any()
        openchromatin_any = sub["openchromatin_name"].fillna("").astype(str).str.strip().ne("").any()

        rows.append(
            {
                "candidate_gene": gene,
                "n_probes": int(len(sub)),
                "best_abs_delta_median": best_abs,
                "min_q_value_bh": min_q,
                "directions": directions,
                "chromosomes": chrs,
                "regulatory_feature_groups": reg_groups_str,
                "tfbs_any": bool(tfbs_any),
                "enhancer_any": bool(enhancer_any),
                "dnase_any": bool(dnase_any),
                "openchromatin_any": bool(openchromatin_any),
            }
        )

    out = pd.DataFrame(rows).sort_values(
        ["n_probes", "best_abs_delta_median", "candidate_gene"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    return out


def build_regulatory_feature_overlap(hits: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    if hits.empty:
        return pd.DataFrame(columns=["candidate_gene", "feature_type", "feature_value", "n_probes"])

    rows = []
    feature_cols = {
        "regulatory_feature_group": "regulatory_feature_group",
        "regulatory_feature_name": "regulatory_feature_name",
        "phantom4_enhancers": "phantom4_enhancers",
        "dnase_name": "dnase_name",
        "openchromatin_name": "openchromatin_name",
        "tfbs_name": "tfbs_name",
    }

    for gene in genes:
        sub = hits[hits["candidate_gene"] == gene].copy()
        if sub.empty:
            continue

        for feature_type, col in feature_cols.items():
            vals = sub[col].fillna("").astype(str).str.strip()
            vals = vals[vals != ""]
            if vals.empty:
                continue
            counts = vals.value_counts()
            for feature_value, n in counts.items():
                rows.append(
                    {
                        "candidate_gene": gene,
                        "feature_type": feature_type,
                        "feature_value": feature_value,
                        "n_probes": int(n),
                    }
                )

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    out = out.sort_values(
        ["feature_type", "n_probes", "candidate_gene", "feature_value"],
        ascending=[True, False, True, True],
    ).reset_index(drop=True)
    return out


def build_top200_intersection(full_hits: pd.DataFrame, top200_hits: pd.DataFrame) -> pd.DataFrame:
    genes = sorted(set(full_hits["candidate_gene"].tolist()))
    rows = []
    for gene in genes:
        full_sub = full_hits[full_hits["candidate_gene"] == gene]
        top_sub = top200_hits[top200_hits["candidate_gene"] == gene] if not top200_hits.empty else pd.DataFrame()

        rows.append(
            {
                "candidate_gene": gene,
                "n_probes_full_universe": int(len(full_sub)),
                "n_probes_top200": int(len(top_sub)),
                "present_in_top200": bool(len(top_sub) > 0),
                "best_abs_delta_median_full_universe": pd.to_numeric(full_sub["abs_delta_median"], errors="coerce").max(),
                "best_abs_delta_median_top200": pd.to_numeric(top_sub["abs_delta_median"], errors="coerce").max() if len(top_sub) > 0 else None,
            }
        )

    out = pd.DataFrame(rows).sort_values(
        ["present_in_top200", "n_probes_full_universe", "candidate_gene"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    return out


def write_candidate_summary_tex(df: pd.DataFrame, outpath: Path) -> None:
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Focused evidence summary for candidate genes detected in the full-universe D3 versus not-D3 methylation scan.}")
    lines.append(r"\label{tab:focused_candidate_evidence}")
    lines.append(r"\begin{tabular}{lccccc}")
    lines.append(r"\hline")
    lines.append(r"Gene & Probes & Best $|\Delta \mathrm{median}|$ & Min $q$ & TFBS & Enhancer \\")
    lines.append(r"\hline")

    for _, row in df.iterrows():
        gene = latex_escape(row["candidate_gene"])
        n_probes = latex_escape(row["n_probes"])
        best_abs = "" if pd.isna(row["best_abs_delta_median"]) else f"{float(row['best_abs_delta_median']):.4f}"
        min_q = "" if pd.isna(row["min_q_value_bh"]) else f"{float(row['min_q_value_bh']):.4f}"
        tfbs_any = "yes" if bool(row["tfbs_any"]) else "no"
        enhancer_any = "yes" if bool(row["enhancer_any"]) else "no"
        lines.append(f"{gene} & {n_probes} & {best_abs} & {min_q} & {tfbs_any} & {enhancer_any} \\\\")

    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    outpath.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_top200_intersection_tex(df: pd.DataFrame, outpath: Path) -> None:
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Intersection of focused candidate genes with the original top-200 directional methylation export.}")
    lines.append(r"\label{tab:focused_candidate_top200_intersection}")
    lines.append(r"\begin{tabular}{lccc}")
    lines.append(r"\hline")
    lines.append(r"Gene & Full-universe probes & Top-200 probes & Present in top-200 \\")
    lines.append(r"\hline")

    for _, row in df.iterrows():
        gene = latex_escape(row["candidate_gene"])
        n_full = latex_escape(row["n_probes_full_universe"])
        n_top = latex_escape(row["n_probes_top200"])
        present = "yes" if bool(row["present_in_top200"]) else "no"
        lines.append(f"{gene} & {n_full} & {n_top} & {present} \\\\")

    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    outpath.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_effect_plot(df: pd.DataFrame, outpath: Path) -> None:
    if df.empty:
        return

    tmp = df.copy()
    tmp = tmp.sort_values(["candidate_gene", "best_abs_delta_median"], ascending=[True, False]).reset_index(drop=True)
    tmp["rank"] = range(1, len(tmp) + 1)

    plt.figure(figsize=(9.5, 4.8))
    plt.scatter(tmp["rank"], tmp["best_abs_delta_median"], alpha=0.9)
    for _, row in tmp.iterrows():
        plt.annotate(
            str(row["candidate_gene"]),
            (row["rank"], row["best_abs_delta_median"]),
            fontsize=8,
            alpha=0.8,
        )
    plt.title("Focused candidate evidence, best absolute methylation effect")
    plt.xlabel("Candidate rank")
    plt.ylabel("Best |delta median|")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def write_logbook_entry(summary: dict, outpath: Path) -> None:
    shared_reg_groups = summary["shared_regulatory_feature_groups"]
    shared_reg_groups_text = ", ".join(shared_reg_groups) if shared_reg_groups else "none detected"

    genes_in_top200 = summary["genes_present_in_top200"]
    genes_in_top200_text = ", ".join(genes_in_top200) if genes_in_top200 else "none"

    tex = rf"""
\subsection*{{Logbook Entry {ENTRY_NUMBER}: Focused candidate evidence and regulatory-overlap follow up}}

\textbf{{What was being measured.}}
A focused follow up was performed on the six candidate genes detected in the full-universe D3 versus not-D3 methylation scan, with the aim of comparing their probe-level evidence, regulatory annotations, and relationship to the original top-200 directional export. The broader question was whether the SIM2-centered 21q signature sits inside a more coherent oncogenic developmental program.

\textbf{{Methods.}}
Probe-level rows were extracted from the full-universe D3 versus not-D3 table for PAX6, GLI2, BMP4, PTCH1, SHH, and OLIG2, after manifest-based annotation merge. For each candidate, the analysis summarized probe count, best absolute median difference, minimum adjusted \(q\)-value, chromosome, genomic context, and available regulatory annotations including regulatory feature group, regulatory feature name, PHANTOM4 enhancer labels, DNase labels, open chromatin labels, and TFBS labels. The same candidate set was then intersected with the original top-200 directional methylation export to determine which signals rose to the shortlist level.

\textbf{{Results.}}
The full-universe candidate set contained probe-level support for PAX6, GLI2, BMP4, PTCH1, SHH, and OLIG2. Shared regulatory feature groups across these genes were: {latex_escape(shared_reg_groups_text)}. Intersection with the original top-200 directional export showed top-200 presence for: {latex_escape(genes_in_top200_text)}. This establishes which components of the broader developmental bridge were strong enough to enter the original shortlist and which remained weaker background candidates visible only in the broader universe scan.

\textbf{{Interpretation.}}
This step helps determine whether the SIM2-centered signal can reasonably be framed as part of a broader oncogenic developmental program. The logic is not merely nominal overlap of gene names, but convergence across methylation effect size, regulatory-context annotation, and shortlist prominence. Shared regulatory annotations such as enhancer, TFBS, DNase, or open-chromatin labels would strengthen the case for coordinated regulatory structure, whereas scattered weak hits without shared regulatory context would argue for a looser developmental backdrop rather than a coherent program.

\textbf{{Tables written.}}
\texttt{{results/gse240704/focused\_candidate\_evidence/focused\_candidate\_full\_universe\_hits.csv}}\\
\texttt{{results/gse240704/focused\_candidate\_evidence/focused\_candidate\_summary.csv}}\\
\texttt{{results/gse240704/focused\_candidate\_evidence/focused\_candidate\_regulatory\_overlap.csv}}\\
\texttt{{results/gse240704/focused\_candidate\_evidence/focused\_candidate\_top200\_intersection.csv}}\\
\texttt{{results/gse240704/focused\_candidate\_evidence/tables\_tex/table\_focused\_candidate\_evidence.tex}}\\
\texttt{{results/gse240704/focused\_candidate\_evidence/tables\_tex/table\_focused\_candidate\_top200\_intersection.tex}}\\
\texttt{{results/gse240704/focused\_candidate\_evidence/focused\_candidate\_summary.json}}

\textbf{{Figures.}}
\texttt{{results/gse240704/focused\_candidate\_evidence/plots/focused\_candidate\_best\_effect.png}}

\textbf{{Next steps.}}
The next step is to use the focused candidate evidence to decide whether the manuscript should describe the broader background as a coordinated developmental program or as a weaker distributed backdrop around the dominant local SIM2 signal. If regulatory-context sharing is sparse, interpretation should remain cautious. If shared enhancer, TFBS, DNase, or shortlist enrichment is clear across several genes, the developmental-program language becomes more defensible.
""".strip() + "\n"

    outpath.write_text(tex, encoding="utf-8")


def main() -> None:
    ensure_dir(OUTDIR)
    ensure_dir(PLOTDIR)
    ensure_dir(TABLETEXDIR)
    ensure_dir(LOGBOOKDIR)

    full_df = load_table_with_manifest(FULL_UNIVERSE_TABLE)
    top_df = load_table_with_manifest(TOP200_DIRECTIONAL_TABLE)

    focused_full_hits = extract_gene_probe_rows(full_df, FOCUSED_GENES)
    focused_top_hits = extract_gene_probe_rows(top_df, FOCUSED_GENES)

    full_hits_path = OUTDIR / "focused_candidate_full_universe_hits.csv"
    focused_full_hits.to_csv(full_hits_path, index=False)

    summary_df = summarize_by_gene(focused_full_hits)
    summary_path = OUTDIR / "focused_candidate_summary.csv"
    summary_df.to_csv(summary_path, index=False)

    reg_overlap_df = build_regulatory_feature_overlap(focused_full_hits, FOCUSED_GENES)
    reg_overlap_path = OUTDIR / "focused_candidate_regulatory_overlap.csv"
    reg_overlap_df.to_csv(reg_overlap_path, index=False)

    top200_intersection_df = build_top200_intersection(focused_full_hits, focused_top_hits)
    top200_intersection_path = OUTDIR / "focused_candidate_top200_intersection.csv"
    top200_intersection_df.to_csv(top200_intersection_path, index=False)

    # Shared regulatory feature groups across >=2 genes
    shared_reg_groups = []
    if not reg_overlap_df.empty:
        rg = reg_overlap_df[reg_overlap_df["feature_type"] == "regulatory_feature_group"].copy()
        if not rg.empty:
            counts = rg.groupby("feature_value")["candidate_gene"].nunique().sort_values(ascending=False)
            shared_reg_groups = counts[counts >= 2].index.astype(str).tolist()

    genes_in_top200 = top200_intersection_df.loc[
        top200_intersection_df["present_in_top200"], "candidate_gene"
    ].tolist()

    table_evidence_tex = TABLETEXDIR / "table_focused_candidate_evidence.tex"
    write_candidate_summary_tex(summary_df, table_evidence_tex)

    table_top200_tex = TABLETEXDIR / "table_focused_candidate_top200_intersection.tex"
    write_top200_intersection_tex(top200_intersection_df, table_top200_tex)

    plot_path = PLOTDIR / "focused_candidate_best_effect.png"
    write_effect_plot(summary_df, plot_path)

    summary_json = {
        "dataset": DATASET,
        "focused_genes": FOCUSED_GENES,
        "full_universe_table": str(FULL_UNIVERSE_TABLE),
        "top200_directional_table": str(TOP200_DIRECTIONAL_TABLE),
        "n_focused_genes": len(FOCUSED_GENES),
        "n_detected_focused_genes": int(len(summary_df)),
        "genes_present_in_top200": genes_in_top200,
        "shared_regulatory_feature_groups": shared_reg_groups,
        "outputs": {
            "full_hits_csv": str(full_hits_path),
            "summary_csv": str(summary_path),
            "regulatory_overlap_csv": str(reg_overlap_path),
            "top200_intersection_csv": str(top200_intersection_path),
            "evidence_tex": str(table_evidence_tex),
            "top200_tex": str(table_top200_tex),
            "plot": str(plot_path),
        },
    }

    summary_json_path = OUTDIR / "focused_candidate_summary.json"
    with open(summary_json_path, "w", encoding="utf-8") as f:
        json.dump(summary_json, f, indent=2)

    logbook_path = LOGBOOKDIR / "entry_31_focused_candidate_evidence.tex"
    write_logbook_entry(summary_json, logbook_path)

    print(json.dumps(summary_json, indent=2))
    print(f"[ok] wrote logbook entry to {logbook_path}")


if __name__ == "__main__":
    main()
