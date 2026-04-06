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

DIRECTIONAL_TABLE = (
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

SIM2_SEED_TABLE = (
    REPO_ROOT
    / "results"
    / DATASET
    / "sim2_followup"
    / "SIM2_downstream_target_seed_table.csv"
)

OUTDIR = (
    REPO_ROOT
    / "results"
    / DATASET
    / "candidate_gene_directional_check"
)

PLOTDIR = OUTDIR / "plots"
TABLETEXDIR = OUTDIR / "tables_tex"
LOGBOOKDIR = OUTDIR / "logbook"

ENTRY_NUMBER = 27

# Core requested set
CORE_GENES = ["PTCH1", "GLI1", "GLI2", "GLI3", "SHH", "SOX2", "OLIG2"]

# Also fold in the rest of the SIM2 seed list if present
EXTRA_GENES = ["BMP4", "PAX6", "ASCL1", "MYCN"]


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


def load_candidate_genes() -> list[str]:
    genes = []
    genes.extend(CORE_GENES)
    genes.extend(EXTRA_GENES)

    if SIM2_SEED_TABLE.exists():
        seed = pd.read_csv(SIM2_SEED_TABLE)
        if "gene" in seed.columns:
            genes.extend(seed["gene"].dropna().astype(str).tolist())

    # unique, preserve order
    seen = set()
    ordered = []
    for g in genes:
        gg = g.strip().upper()
        if gg and gg not in seen:
            ordered.append(gg)
            seen.add(gg)
    return ordered


def build_merged_table() -> pd.DataFrame:
    if not DIRECTIONAL_TABLE.exists():
        raise FileNotFoundError(f"Missing directional table: {DIRECTIONAL_TABLE}")
    if not MANIFEST_TABLE.exists():
        raise FileNotFoundError(f"Missing manifest table: {MANIFEST_TABLE}")

    directional = normalize_columns(pd.read_csv(DIRECTIONAL_TABLE))
    manifest = normalize_columns(pd.read_parquet(MANIFEST_TABLE))

    directional = deduplicate_columns(directional)
    manifest = deduplicate_columns(manifest)

    probe_col_dir = find_first_existing(
        directional, ["probe_id", "ilmnid", "ilmn_id", "name", "id_ref", "id"]
    )
    probe_col_ann = find_first_existing(
        manifest, ["probe_id", "ilmnid", "ilmn_id", "name", "id_ref", "id"]
    )
    if probe_col_dir is None or probe_col_ann is None:
        raise ValueError("Could not identify probe id columns for merge.")

    directional[probe_col_dir] = directional[probe_col_dir].astype(str).str.strip()
    manifest[probe_col_ann] = manifest[probe_col_ann].astype(str).str.strip()

    merged = directional.merge(
        manifest,
        left_on=probe_col_dir,
        right_on=probe_col_ann,
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


def choose_effect_col(df: pd.DataFrame) -> Optional[str]:
    return find_first_existing(
        df,
        [
            "delta_median_a_minus_b",
            "delta_mean_a_minus_b",
            "abs_delta_median",
            "rank_biserial",
            "abs_rank_biserial",
        ],
    )


def choose_direction_col(df: pd.DataFrame) -> Optional[str]:
    return find_first_existing(df, ["direction", "probe_direction", "contrast_direction"])


def choose_pos_col(df: pd.DataFrame) -> Optional[str]:
    return find_first_existing(df, ["mapinfo_ann", "mapinfo_dir", "mapinfo", "position_std", "position"])


def choose_chr_col(df: pd.DataFrame) -> Optional[str]:
    return find_first_existing(df, ["chromosome_ann", "chromosome_dir", "chromosome", "chromosome_std", "chr"])


def choose_context_cols(df: pd.DataFrame) -> tuple[Optional[str], Optional[str]]:
    refgene_group = find_first_existing(
        df,
        [
            "refgene_group_ann",
            "refgene_group_dir",
            "ucsc_refgene_group_ann",
            "ucsc_refgene_group_dir",
            "refgene_group",
        ],
    )
    island = find_first_existing(
        df,
        [
            "relation_to_cpg_island_ann",
            "relation_to_cpg_island_dir",
            "relation_to_ucsc_cpg_island_ann",
            "relation_to_ucsc_cpg_island_dir",
            "relation_to_cpg_island",
        ],
    )
    return refgene_group, island


def build_candidate_hit_table(merged: pd.DataFrame, candidate_genes: list[str]) -> pd.DataFrame:
    gene_col = choose_gene_col(merged)
    if gene_col is None:
        raise ValueError("Could not identify gene annotation column.")

    effect_col = choose_effect_col(merged)
    direction_col = choose_direction_col(merged)
    pos_col = choose_pos_col(merged)
    chr_col = choose_chr_col(merged)
    refgene_group_col, island_col = choose_context_cols(merged)

    rows = []
    for _, row in merged.iterrows():
        tokens = {tok.upper() for tok in split_gene_tokens(row.get(gene_col))}
        overlap = sorted(tokens.intersection(candidate_genes))
        if not overlap:
            continue

        for g in overlap:
            rows.append(
                {
                    "candidate_gene": g,
                    "probe_id": row.get("probe_id"),
                    "direction": row.get(direction_col) if direction_col else None,
                    "delta_median_a_minus_b": row.get("delta_median_a_minus_b"),
                    "delta_mean_a_minus_b": row.get("delta_mean_a_minus_b"),
                    "abs_delta_median": row.get("abs_delta_median"),
                    "rank_biserial": row.get("rank_biserial"),
                    "q_value_bh": row.get("q_value_bh"),
                    "p_value": row.get("p_value"),
                    "chromosome": row.get(chr_col) if chr_col else None,
                    "position": row.get(pos_col) if pos_col else None,
                    "refgene_group": row.get(refgene_group_col) if refgene_group_col else None,
                    "relation_to_cpg_island": row.get(island_col) if island_col else None,
                    "raw_gene_annotation": row.get(gene_col),
                }
            )

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    if "abs_delta_median" in out.columns:
        out["abs_delta_median_num"] = pd.to_numeric(out["abs_delta_median"], errors="coerce")
        out = out.sort_values(
            ["candidate_gene", "abs_delta_median_num"],
            ascending=[True, False],
        ).drop(columns=["abs_delta_median_num"])
    else:
        out = out.sort_values(["candidate_gene", "probe_id"])

    return out.reset_index(drop=True)


def build_candidate_summary(hits: pd.DataFrame, candidate_genes: list[str]) -> pd.DataFrame:
    rows = []
    genes_present = set(hits["candidate_gene"].unique()) if not hits.empty else set()

    for g in candidate_genes:
        if g not in genes_present:
            rows.append(
                {
                    "candidate_gene": g,
                    "n_probes": 0,
                    "directions": "",
                    "best_abs_delta_median": None,
                    "min_q_value_bh": None,
                    "status": "not_detected_in_directional_table",
                }
            )
            continue

        sub = hits[hits["candidate_gene"] == g].copy()
        best_abs = pd.to_numeric(sub["abs_delta_median"], errors="coerce").max()
        min_q = pd.to_numeric(sub["q_value_bh"], errors="coerce").min()
        directions = ",".join(sorted(set(sub["direction"].dropna().astype(str).tolist())))
        rows.append(
            {
                "candidate_gene": g,
                "n_probes": int(len(sub)),
                "directions": directions,
                "best_abs_delta_median": best_abs,
                "min_q_value_bh": min_q,
                "status": "detected_in_directional_table",
            }
        )

    out = pd.DataFrame(rows)
    out = out.sort_values(
        ["status", "n_probes", "best_abs_delta_median", "candidate_gene"],
        ascending=[True, False, False, True],
    ).reset_index(drop=True)
    return out


def write_summary_tex(df: pd.DataFrame, outpath: Path) -> None:
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Directional methylation check for candidate genes linked to the SIM2-centered interpretation, including Hedgehog-related and stemness-related genes. Detection here refers to presence in the exported D3 versus not-D3 directional probe table.}")
    lines.append(r"\label{tab:sim2_candidate_directional_check}")
    lines.append(r"\begin{tabular}{lclll}")
    lines.append(r"\hline")
    lines.append(r"Gene & Probes & Direction & Best $|\Delta \mathrm{median}|$ & Status \\")
    lines.append(r"\hline")

    for _, row in df.iterrows():
        gene = latex_escape(row["candidate_gene"])
        n_probes = latex_escape(row["n_probes"])
        directions = latex_escape(row["directions"])
        best_abs = "" if pd.isna(row["best_abs_delta_median"]) else f"{float(row['best_abs_delta_median']):.4f}"
        best_abs = latex_escape(best_abs)
        status = latex_escape(row["status"])
        lines.append(f"{gene} & {n_probes} & {directions} & {best_abs} & {status} \\\\")

    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    outpath.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_plot(summary: pd.DataFrame, outpath: Path) -> None:
    tmp = summary.copy()
    tmp["n_probes_num"] = pd.to_numeric(tmp["n_probes"], errors="coerce").fillna(0)
    tmp = tmp.sort_values(["n_probes_num", "candidate_gene"], ascending=[False, True]).reset_index(drop=True)
    tmp["rank"] = range(1, len(tmp) + 1)

    plt.figure(figsize=(10.5, 4.8))
    plt.scatter(tmp["rank"], tmp["n_probes_num"], alpha=0.9)
    for _, row in tmp.iterrows():
        plt.annotate(
            str(row["candidate_gene"]),
            (row["rank"], row["n_probes_num"]),
            fontsize=8,
            alpha=0.8,
        )
    plt.title("Candidate-gene directional methylation support in D3 versus not-D3")
    plt.xlabel("Candidate gene rank")
    plt.ylabel("Number of detected probes")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def write_logbook_entry(summary: dict, outpath: Path) -> None:
    detected = summary["detected_genes"]
    not_detected = summary["not_detected_genes"]

    detected_text = ", ".join(detected) if detected else "none"
    not_detected_text = ", ".join(not_detected) if not_detected else "none"

    tex = rf"""
\subsection*{{Logbook Entry {ENTRY_NUMBER}: Directional methylation check for Hedgehog- and stemness-linked candidate genes}}

\textbf{{What was being measured.}}
A candidate-gene overlay was performed to test whether genes related to the SIM2-centered interpretation, especially Hedgehog-linked and stemness-linked genes, were actually represented in the exported D3 versus not-D3 directional methylation signal.

\textbf{{Methods.}}
The directional probe table
\texttt{{results/gse240704/directional\_biology\_followup/D3\_vs\_not\_D3\_top\_200\_probes\_with\_direction.csv}}
was merged with the manifest-derived probe annotation table
\texttt{{results/gse240704/manual\_manifest\_annotation/gpl23976\_annotation\_from\_manifest.parquet}}.
Candidate genes included PTCH1, GLI1, GLI2, GLI3, SHH, SOX2, OLIG2, and additional genes from the curated SIM2 seed list. Probe-level overlaps were identified by manifest gene-token matching. A per-gene summary table, a probe-level hit table, a LaTeX table, and a support plot were exported.

\textbf{{Results.}}
Within the exported directional methylation table, the following candidate genes were detected: {latex_escape(detected_text)}.
The following candidates were not detected in the same directional table: {latex_escape(not_detected_text)}.
This result provides an initial check of whether the SIM2-centered interpretation is accompanied by explicit methylation support for Hedgehog-related or stemness-related genes in the current directional export layer. Because this analysis is restricted to the already exported directional table, absence here should not yet be interpreted as genome-wide absence in the full dataset.

\textbf{{Tables written.}}
\texttt{{results/gse240704/candidate\_gene\_directional\_check/candidate\_gene\_directional\_hits.csv}}\\
\texttt{{results/gse240704/candidate\_gene\_directional\_check/candidate\_gene\_directional\_summary.csv}}\\
\texttt{{results/gse240704/candidate\_gene\_directional\_check/tables\_tex/table\_sim2\_candidate\_directional\_check.tex}}\\
\texttt{{results/gse240704/candidate\_gene\_directional\_check/candidate\_gene\_directional\_summary.json}}

\textbf{{Figures.}}
\texttt{{results/gse240704/candidate\_gene\_directional\_check/plots/candidate\_gene\_directional\_support.png}}

\textbf{{Interpretation.}}
This step provides the first bridge between the SIM2-centered 21q methylation neighborhood and the broader developmental or stemness-oriented interpretation layer. Positive detection of candidates such as PTCH1, GLI-family genes, SHH, SOX2, or OLIG2 in the directional probe export would support continuity between the local SIM2 signature and broader developmental programs. Non-detection in this exported directional set does not eliminate such programs from the dataset, but indicates that they are not among the currently exported directional methylation highlights.

\textbf{{Caveats.}}
This is a candidate-gene overlay on the exported D3 directional table, not a full probe-universe test across all measured methylation sites. A fuller absence or presence statement would require re-querying the complete D3 versus not-D3 differential table if available.

\textbf{{Next steps.}}
If any Hedgehog-linked or stemness-linked candidates are detected, the next step is to compare their directional behavior and effect sizes against the SIM2-centered signal. If none are detected, the next step is to treat the current evidence as SIM2-localized but not yet extended into a broader developmental methylation module, then decide whether a full-universe candidate scan should be run from the complete differential output.
""".strip() + "\n"

    outpath.write_text(tex, encoding="utf-8")


def main() -> None:
    ensure_dir(OUTDIR)
    ensure_dir(PLOTDIR)
    ensure_dir(TABLETEXDIR)
    ensure_dir(LOGBOOKDIR)

    candidate_genes = load_candidate_genes()
    merged = build_merged_table()

    hits = build_candidate_hit_table(merged, candidate_genes=candidate_genes)
    hits_path = OUTDIR / "candidate_gene_directional_hits.csv"
    hits.to_csv(hits_path, index=False)

    summary = build_candidate_summary(hits, candidate_genes=candidate_genes)
    summary_path_csv = OUTDIR / "candidate_gene_directional_summary.csv"
    summary.to_csv(summary_path_csv, index=False)

    tex_path = TABLETEXDIR / "table_sim2_candidate_directional_check.tex"
    write_summary_tex(summary, tex_path)

    plot_path = PLOTDIR / "candidate_gene_directional_support.png"
    write_plot(summary, plot_path)

    detected_genes = summary.loc[
        summary["status"] == "detected_in_directional_table", "candidate_gene"
    ].tolist()
    not_detected_genes = summary.loc[
        summary["status"] == "not_detected_in_directional_table", "candidate_gene"
    ].tolist()

    summary_json = {
        "dataset": DATASET,
        "candidate_genes_tested": candidate_genes,
        "n_candidate_genes_tested": len(candidate_genes),
        "n_detected_genes": len(detected_genes),
        "detected_genes": detected_genes,
        "not_detected_genes": not_detected_genes,
        "outputs": {
            "hits_csv": str(hits_path),
            "summary_csv": str(summary_path_csv),
            "summary_tex": str(tex_path),
            "plot": str(plot_path),
        },
    }

    summary_path_json = OUTDIR / "candidate_gene_directional_summary.json"
    with open(summary_path_json, "w", encoding="utf-8") as f:
        json.dump(summary_json, f, indent=2)

    logbook_path = LOGBOOKDIR / "entry_27_candidate_gene_directional_check.tex"
    write_logbook_entry(summary_json, logbook_path)

    print(json.dumps(summary_json, indent=2))
    print(f"[ok] wrote logbook entry to {logbook_path}")


if __name__ == "__main__":
    main()
