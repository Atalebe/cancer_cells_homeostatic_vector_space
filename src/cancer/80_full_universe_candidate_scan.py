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
    / "full_universe_candidate_scan"
)

PLOTDIR = OUTDIR / "plots"
TABLETEXDIR = OUTDIR / "tables_tex"
LOGBOOKDIR = OUTDIR / "logbook"

ENTRY_NUMBER = 29

CORE_GENES = ["PTCH1", "GLI1", "GLI2", "GLI3", "SHH", "SOX2", "OLIG2"]
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

    seen = set()
    ordered = []
    for g in genes:
        gg = g.strip().upper()
        if gg and gg not in seen:
            ordered.append(gg)
            seen.add(gg)
    return ordered


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


def choose_effect_cols(df: pd.DataFrame) -> list[str]:
    candidates = [
        "delta_median_a_minus_b",
        "delta_mean_a_minus_b",
        "abs_delta_median",
        "rank_biserial",
        "q_value_bh",
        "p_value",
        "u_statistic",
        "mean_a",
        "mean_b",
        "median_a",
        "median_b",
        "std_a",
        "std_b",
    ]
    return [c for c in candidates if c in df.columns]


def choose_direction_col(df: pd.DataFrame) -> Optional[str]:
    return find_first_existing(df, ["direction", "probe_direction", "contrast_direction"])


def choose_context_cols(df: pd.DataFrame) -> tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
    chr_col = find_first_existing(df, ["chromosome_ann", "chromosome_dir", "chromosome", "chromosome_std", "chr"])
    pos_col = find_first_existing(df, ["mapinfo_ann", "mapinfo_dir", "mapinfo", "position_std", "position"])
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
    return chr_col, pos_col, refgene_group, island


def find_universe_table() -> Path:
    candidates = [
        REPO_ROOT / "results" / DATASET / "directional_biology_followup" / "D3_vs_not_D3_all_probes_with_direction.csv",
        REPO_ROOT / "results" / DATASET / "directional_biology_followup" / "D3_vs_not_D3_all_probes.csv",
        REPO_ROOT / "results" / DATASET / "directional_biology_followup" / "D3_vs_not_D3_probe_table.csv",
        REPO_ROOT / "results" / DATASET / "d3_interpretation_followup" / "D3_vs_not_D3_all_probes.csv",
        REPO_ROOT / "results" / DATASET / "d3_interpretation_followup" / "D3_vs_not_D3_probe_signatures.csv",
        REPO_ROOT / "results" / DATASET / "d2_d3_probe_signatures" / "D3_vs_not_D3_all_probes.csv",
        REPO_ROOT / "results" / DATASET / "d2_d3_probe_signatures" / "D3_vs_not_D3_probe_signatures.csv",
        REPO_ROOT / "results" / DATASET / "manifest_reannotated" / "D3_vs_not_D3_top_200_probes_normalized_annotation.csv",
        REPO_ROOT / "results" / DATASET / "directional_biology_followup" / "D3_vs_not_D3_top_200_probes_with_direction.csv",
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(
        "Could not locate a plausible D3 vs not-D3 universe table. "
        "Edit find_universe_table() with the correct repo path."
    )


def build_merged_table(universe_table: Path) -> pd.DataFrame:
    if not universe_table.exists():
        raise FileNotFoundError(f"Missing universe table: {universe_table}")
    if not MANIFEST_TABLE.exists():
        raise FileNotFoundError(f"Missing manifest table: {MANIFEST_TABLE}")

    universe = normalize_columns(pd.read_csv(universe_table))
    manifest = normalize_columns(pd.read_parquet(MANIFEST_TABLE))

    universe = deduplicate_columns(universe)
    manifest = deduplicate_columns(manifest)

    probe_col_universe = find_first_existing(
        universe, ["probe_id", "ilmnid", "ilmn_id", "name", "id_ref", "id"]
    )
    probe_col_ann = find_first_existing(
        manifest, ["probe_id", "ilmnid", "ilmn_id", "name", "id_ref", "id"]
    )

    if probe_col_universe is None or probe_col_ann is None:
        raise ValueError("Could not identify probe id columns for merge.")

    universe[probe_col_universe] = universe[probe_col_universe].astype(str).str.strip()
    manifest[probe_col_ann] = manifest[probe_col_ann].astype(str).str.strip()

    merged = universe.merge(
        manifest,
        left_on=probe_col_universe,
        right_on=probe_col_ann,
        how="left",
        suffixes=("_dir", "_ann"),
    )
    merged = normalize_columns(merged)
    merged = deduplicate_columns(merged)
    return merged


def build_candidate_hits(merged: pd.DataFrame, candidate_genes: list[str]) -> pd.DataFrame:
    gene_col = choose_gene_col(merged)
    if gene_col is None:
        raise ValueError("Could not identify gene annotation column.")

    effect_cols = choose_effect_cols(merged)
    direction_col = choose_direction_col(merged)
    chr_col, pos_col, refgene_group_col, island_col = choose_context_cols(merged)

    rows = []
    for _, row in merged.iterrows():
        tokens = {tok.upper() for tok in split_gene_tokens(row.get(gene_col))}
        overlap = sorted(tokens.intersection(candidate_genes))
        if not overlap:
            continue

        for g in overlap:
            out_row = {
                "candidate_gene": g,
                "probe_id": row.get("probe_id"),
                "direction": row.get(direction_col) if direction_col else None,
                "chromosome": row.get(chr_col) if chr_col else None,
                "position": row.get(pos_col) if pos_col else None,
                "refgene_group": row.get(refgene_group_col) if refgene_group_col else None,
                "relation_to_cpg_island": row.get(island_col) if island_col else None,
                "raw_gene_annotation": row.get(gene_col),
            }
            for c in effect_cols:
                out_row[c] = row.get(c)
            rows.append(out_row)

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    sort_col = "abs_delta_median" if "abs_delta_median" in out.columns else None
    if sort_col is None and "delta_median_a_minus_b" in out.columns:
        out["_abs_sort"] = pd.to_numeric(out["delta_median_a_minus_b"], errors="coerce").abs()
        out = out.sort_values(["candidate_gene", "_abs_sort"], ascending=[True, False]).drop(columns=["_abs_sort"])
    elif sort_col is not None:
        out["_abs_sort"] = pd.to_numeric(out[sort_col], errors="coerce")
        out = out.sort_values(["candidate_gene", "_abs_sort"], ascending=[True, False]).drop(columns=["_abs_sort"])
    else:
        out = out.sort_values(["candidate_gene", "probe_id"], ascending=[True, True])

    return out.reset_index(drop=True)


def build_candidate_summary(hits: pd.DataFrame, candidate_genes: list[str]) -> pd.DataFrame:
    rows = []
    detected = set(hits["candidate_gene"].unique()) if not hits.empty else set()

    for g in candidate_genes:
        if g not in detected:
            rows.append(
                {
                    "candidate_gene": g,
                    "n_probes": 0,
                    "directions": "",
                    "best_abs_delta_median": None,
                    "min_q_value_bh": None,
                    "chromosomes": "",
                    "status": "not_detected_in_universe_table",
                }
            )
            continue

        sub = hits[hits["candidate_gene"] == g].copy()
        best_abs = None
        if "abs_delta_median" in sub.columns:
            best_abs = pd.to_numeric(sub["abs_delta_median"], errors="coerce").max()
        elif "delta_median_a_minus_b" in sub.columns:
            best_abs = pd.to_numeric(sub["delta_median_a_minus_b"], errors="coerce").abs().max()

        min_q = None
        if "q_value_bh" in sub.columns:
            min_q = pd.to_numeric(sub["q_value_bh"], errors="coerce").min()

        directions = ",".join(sorted(set(sub["direction"].dropna().astype(str).tolist()))) if "direction" in sub.columns else ""
        chrs = ",".join(sorted(set(sub["chromosome"].dropna().astype(str).tolist()))) if "chromosome" in sub.columns else ""

        rows.append(
            {
                "candidate_gene": g,
                "n_probes": int(len(sub)),
                "directions": directions,
                "best_abs_delta_median": best_abs,
                "min_q_value_bh": min_q,
                "chromosomes": chrs,
                "status": "detected_in_universe_table",
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
    lines.append(r"\caption{Full-universe candidate scan for Hedgehog-linked, stemness-linked, and SIM2-seed genes in the D3 versus not-D3 methylation comparison. Detection here refers to presence in the selected D3 versus not-D3 probe universe table used for this scan.}")
    lines.append(r"\label{tab:full_universe_candidate_scan}")
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


def write_support_plot(summary: pd.DataFrame, outpath: Path) -> None:
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
    plt.title("Full-universe candidate-gene support in D3 versus not-D3")
    plt.xlabel("Candidate gene rank")
    plt.ylabel("Number of detected probes")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def write_logbook_entry(summary_json: dict, outpath: Path) -> None:
    detected = ", ".join(summary_json["detected_genes"]) if summary_json["detected_genes"] else "none"
    not_detected = ", ".join(summary_json["not_detected_genes"]) if summary_json["not_detected_genes"] else "none"
    universe_path = summary_json["universe_table_used"]

    tex = rf"""
\subsection*{{Logbook Entry {ENTRY_NUMBER}: Full-universe candidate scan for Hedgehog- and stemness-linked genes}}

\textbf{{What was being measured.}}
A full-universe candidate scan was performed to determine whether Hedgehog-linked, stemness-linked, and SIM2-seed candidate genes were represented anywhere in the selected D3 versus not-D3 methylation probe universe, rather than only in the previously exported directional highlight table.

\textbf{{Methods.}}
A D3 versus not-D3 probe universe table was identified at
\texttt{{{latex_escape(universe_path)}}}
and merged with the manifest-derived annotation table
\texttt{{results/gse240704/manual\_manifest\_annotation/gpl23976\_annotation\_from\_manifest.parquet}}.
Candidate genes included PTCH1, GLI1, GLI2, GLI3, SHH, SOX2, OLIG2, BMP4, PAX6, ASCL1, MYCN, and genes from the SIM2 seed list. Probe-level overlaps were identified by gene-token matching, then summarized by number of detected probes, direction labels when present, and best available absolute effect size.

\textbf{{Results.}}
Detected candidate genes in the selected universe table were: {latex_escape(detected)}.
Candidate genes not detected in the same universe table were: {latex_escape(not_detected)}.
This step extends the earlier candidate overlay beyond the exported directional shortlist and therefore provides a stronger basis for interpreting whether the SIM2-centered methylation signal sits inside a broader developmental, Hedgehog-related, or stemness-related methylation pattern.

\textbf{{Tables written.}}
\texttt{{results/gse240704/full\_universe\_candidate\_scan/full\_universe\_candidate\_hits.csv}}\\
\texttt{{results/gse240704/full\_universe\_candidate\_scan/full\_universe\_candidate\_summary.csv}}\\
\texttt{{results/gse240704/full\_universe\_candidate\_scan/tables\_tex/table\_full\_universe\_candidate\_scan.tex}}\\
\texttt{{results/gse240704/full\_universe\_candidate\_scan/full\_universe\_candidate\_summary.json}}

\textbf{{Figures.}}
\texttt{{results/gse240704/full\_universe\_candidate\_scan/plots/full\_universe\_candidate\_support.png}}

\textbf{{Interpretation.}}
This step determines whether the earlier absence of Hedgehog- and stemness-linked genes from the exported directional highlight set reflects a true broader absence in the selected D3 versus not-D3 universe, or only a shortlist effect. The interpretation should depend on the detected-versus-not-detected pattern in this broader scan, not on the earlier highlight table alone.

\textbf{{Caveats.}}
This result still depends on the specific D3 versus not-D3 universe table used for the scan. If the selected table is not the complete probe-universe output, then absence statements remain conditional on that source file.

\textbf{{Next steps.}}
Use the full-universe scan result to decide whether the manuscript should describe the SIM2-centered signal as largely local and isolated, or whether it sits within a broader candidate-gene methylation program. If the universe table used here is still only partial, the next step is to point this same script to the true complete D3 versus not-D3 differential output.
""".strip() + "\n"

    outpath.write_text(tex, encoding="utf-8")


def main() -> None:
    ensure_dir(OUTDIR)
    ensure_dir(PLOTDIR)
    ensure_dir(TABLETEXDIR)
    ensure_dir(LOGBOOKDIR)

    candidate_genes = load_candidate_genes()
    universe_table = find_universe_table()
    merged = build_merged_table(universe_table)

    hits = build_candidate_hits(merged, candidate_genes)
    hits_path = OUTDIR / "full_universe_candidate_hits.csv"
    hits.to_csv(hits_path, index=False)

    summary = build_candidate_summary(hits, candidate_genes)
    summary_csv = OUTDIR / "full_universe_candidate_summary.csv"
    summary.to_csv(summary_csv, index=False)

    tex_path = TABLETEXDIR / "table_full_universe_candidate_scan.tex"
    write_summary_tex(summary, tex_path)

    plot_path = PLOTDIR / "full_universe_candidate_support.png"
    write_support_plot(summary, plot_path)

    detected = summary.loc[
        summary["status"] == "detected_in_universe_table", "candidate_gene"
    ].tolist()
    not_detected = summary.loc[
        summary["status"] == "not_detected_in_universe_table", "candidate_gene"
    ].tolist()

    summary_json = {
        "dataset": DATASET,
        "universe_table_used": str(universe_table),
        "candidate_genes_tested": candidate_genes,
        "n_candidate_genes_tested": len(candidate_genes),
        "n_detected_genes": len(detected),
        "detected_genes": detected,
        "not_detected_genes": not_detected,
        "outputs": {
            "hits_csv": str(hits_path),
            "summary_csv": str(summary_csv),
            "summary_tex": str(tex_path),
            "plot": str(plot_path),
        },
    }

    summary_json_path = OUTDIR / "full_universe_candidate_summary.json"
    with open(summary_json_path, "w", encoding="utf-8") as f:
        json.dump(summary_json, f, indent=2)

    logbook_path = LOGBOOKDIR / "entry_29_full_universe_candidate_scan.tex"
    write_logbook_entry(summary_json, logbook_path)

    print(json.dumps(summary_json, indent=2))
    print(f"[ok] wrote logbook entry to {logbook_path}")


if __name__ == "__main__":
    main()
