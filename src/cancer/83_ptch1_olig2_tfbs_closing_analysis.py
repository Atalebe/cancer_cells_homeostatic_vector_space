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

FOCUSED_HITS_TABLE = (
    REPO_ROOT
    / "results"
    / DATASET
    / "focused_candidate_evidence"
    / "focused_candidate_full_universe_hits.csv"
)

FOCUSED_SUMMARY_TABLE = (
    REPO_ROOT
    / "results"
    / DATASET
    / "focused_candidate_evidence"
    / "focused_candidate_summary.csv"
)

TOP200_INTERSECTION_TABLE = (
    REPO_ROOT
    / "results"
    / DATASET
    / "focused_candidate_evidence"
    / "focused_candidate_top200_intersection.csv"
)

MATCHED_BACKGROUND_TABLE = (
    REPO_ROOT
    / "results"
    / DATASET
    / "matched_background_regulatory_enrichment"
    / "matched_background_enrichment_summary.csv"
)

OUTDIR = (
    REPO_ROOT
    / "results"
    / DATASET
    / "ptch1_olig2_tfbs_closing_analysis"
)

PLOTDIR = OUTDIR / "plots"
TABLETEXDIR = OUTDIR / "tables_tex"
LOGBOOKDIR = OUTDIR / "logbook"

ENTRY_NUMBER = 35
TARGET_GENES = ["PTCH1", "OLIG2"]


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


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


def looks_like_coordinate(text: str) -> bool:
    s = str(text).strip()
    if not s:
        return False
    return bool(re.fullmatch(r"chr[0-9XYM]+:\d+-\d+", s))


def extract_named_tf_candidates(series: pd.Series) -> list[str]:
    """
    Return non-coordinate TFBS-like strings that could plausibly be named TFs.
    """
    vals = []
    for x in series.fillna("").astype(str):
        s = x.strip()
        if not s:
            continue
        if looks_like_coordinate(s):
            continue
        vals.append(s)
    return sorted(set(vals))


def build_gene_detail_table(df: pd.DataFrame, gene: str) -> pd.DataFrame:
    sub = df[df["candidate_gene"].astype(str).str.upper() == gene.upper()].copy()
    if sub.empty:
        return sub
    if "abs_delta_median" in sub.columns:
        sub["_sort_abs"] = pd.to_numeric(sub["abs_delta_median"], errors="coerce")
        sub = sub.sort_values("_sort_abs", ascending=False).drop(columns=["_sort_abs"])
    return sub.reset_index(drop=True)


def summarize_tfbs_identity(df: pd.DataFrame, gene: str) -> dict:
    sub = build_gene_detail_table(df, gene)
    if sub.empty:
        return {
            "candidate_gene": gene,
            "n_probes": 0,
            "n_tfbs_annotated_probes": 0,
            "tfbs_entries": [],
            "all_tfbs_entries_are_coordinates": None,
            "named_tf_candidates": [],
            "master_tf_inference_status": "gene_not_detected",
        }

    tfbs_series = sub["tfbs_name"] if "tfbs_name" in sub.columns else pd.Series([], dtype=object)
    tfbs_nonempty = tfbs_series.fillna("").astype(str).str.strip()
    tfbs_nonempty = tfbs_nonempty[tfbs_nonempty != ""]

    tfbs_entries = sorted(set(tfbs_nonempty.tolist()))
    named_tf_candidates = extract_named_tf_candidates(tfbs_nonempty)

    all_coords = None
    if len(tfbs_entries) > 0:
        all_coords = all(looks_like_coordinate(x) for x in tfbs_entries)

    if len(tfbs_entries) == 0:
        status = "no_tfbs_annotation_present"
    elif all_coords and len(named_tf_candidates) == 0:
        status = "tfbs_intervals_present_but_no_named_tf_identity"
    elif len(named_tf_candidates) > 0:
        status = "named_tf_candidates_present"
    else:
        status = "tfbs_annotation_ambiguous"

    return {
        "candidate_gene": gene,
        "n_probes": int(len(sub)),
        "n_tfbs_annotated_probes": int(len(tfbs_nonempty)),
        "tfbs_entries": tfbs_entries,
        "all_tfbs_entries_are_coordinates": all_coords,
        "named_tf_candidates": named_tf_candidates,
        "master_tf_inference_status": status,
    }


def build_combined_summary(detail_df: pd.DataFrame, top200_df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for gene in TARGET_GENES:
        sub = build_gene_detail_table(detail_df, gene)
        top_sub = top200_df[top200_df["candidate_gene"].astype(str).str.upper() == gene.upper()].copy()

        tfbs_nonempty = (
            sub["tfbs_name"].fillna("").astype(str).str.strip().replace("", pd.NA).dropna()
            if "tfbs_name" in sub.columns else pd.Series([], dtype=object)
        )
        dnase_nonempty = (
            sub["dnase_name"].fillna("").astype(str).str.strip().replace("", pd.NA).dropna()
            if "dnase_name" in sub.columns else pd.Series([], dtype=object)
        )
        openchrom_nonempty = (
            sub["openchromatin_name"].fillna("").astype(str).str.strip().replace("", pd.NA).dropna()
            if "openchromatin_name" in sub.columns else pd.Series([], dtype=object)
        )
        reg_group_nonempty = (
            sub["regulatory_feature_group"].fillna("").astype(str).str.strip().replace("", pd.NA).dropna()
            if "regulatory_feature_group" in sub.columns else pd.Series([], dtype=object)
        )

        best_abs = pd.to_numeric(sub["abs_delta_median"], errors="coerce").max() if len(sub) else None
        min_q = pd.to_numeric(sub["q_value_bh"], errors="coerce").min() if len(sub) else None

        tfbs_entries = sorted(set(tfbs_nonempty.tolist()))
        named_tf_candidates = extract_named_tf_candidates(tfbs_nonempty)
        all_coords = all(looks_like_coordinate(x) for x in tfbs_entries) if tfbs_entries else None

        rows.append(
            {
                "candidate_gene": gene,
                "n_probes": int(len(sub)),
                "best_abs_delta_median": best_abs,
                "min_q_value_bh": min_q,
                "n_tfbs_annotated_probes": int(len(tfbs_nonempty)),
                "tfbs_entries": "; ".join(tfbs_entries),
                "named_tf_candidates": "; ".join(named_tf_candidates),
                "all_tfbs_entries_are_coordinates": all_coords,
                "n_dnase_annotated_probes": int(len(dnase_nonempty)),
                "n_openchromatin_annotated_probes": int(len(openchrom_nonempty)),
                "regulatory_feature_groups": "; ".join(sorted(set(reg_group_nonempty.tolist()))),
                "present_in_top200": bool(len(top_sub) > 0 and top_sub["present_in_top200"].astype(bool).any()),
                "n_top200_probes": int(top_sub["n_probes_top200"].iloc[0]) if len(top_sub) else 0,
            }
        )

    return pd.DataFrame(rows)


def write_summary_tex(df: pd.DataFrame, outpath: Path) -> None:
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Closing analysis of PTCH1 and OLIG2 regulatory annotation support in the D3 versus not-D3 methylation universe.}")
    lines.append(r"\label{tab:ptch1_olig2_tfbs_closing}")
    lines.append(r"\begin{tabular}{lcccccc}")
    lines.append(r"\hline")
    lines.append(r"Gene & Probes & Best $|\Delta \mathrm{median}|$ & TFBS probes & DNase probes & Open chrom. & Top-200 \\")
    lines.append(r"\hline")

    for _, row in df.iterrows():
        gene = latex_escape(row["candidate_gene"])
        n_probes = latex_escape(row["n_probes"])
        best_abs = "" if pd.isna(row["best_abs_delta_median"]) else f"{float(row['best_abs_delta_median']):.4f}"
        n_tfbs = latex_escape(row["n_tfbs_annotated_probes"])
        n_dnase = latex_escape(row["n_dnase_annotated_probes"])
        n_open = latex_escape(row["n_openchromatin_annotated_probes"])
        top200 = "yes" if bool(row["present_in_top200"]) else "no"
        lines.append(f"{gene} & {n_probes} & {best_abs} & {n_tfbs} & {n_dnase} & {n_open} & {top200} \\\\")

    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    outpath.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_effect_plot(df: pd.DataFrame, outpath: Path) -> None:
    if df.empty:
        return

    tmp = df.copy()
    tmp["best_abs_delta_median"] = pd.to_numeric(tmp["best_abs_delta_median"], errors="coerce")
    tmp = tmp.sort_values("best_abs_delta_median", ascending=False).reset_index(drop=True)
    tmp["rank"] = range(1, len(tmp) + 1)

    plt.figure(figsize=(7.5, 4.5))
    plt.scatter(tmp["rank"], tmp["best_abs_delta_median"], alpha=0.9)
    for _, row in tmp.iterrows():
        plt.annotate(
            str(row["candidate_gene"]),
            (row["rank"], row["best_abs_delta_median"]),
            fontsize=8,
            alpha=0.8,
        )
    plt.title("PTCH1 and OLIG2, best absolute methylation effect")
    plt.xlabel("Gene rank")
    plt.ylabel("Best |delta median|")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def write_logbook_entry(summary_json: dict, outpath: Path) -> None:
    master_status = summary_json["master_tf_inference_overall"]
    tfbs_note = summary_json["tfbs_identity_statement"]

    tex = rf"""
\subsection*{{Logbook Entry {ENTRY_NUMBER}: Closing PTCH1 and OLIG2 TFBS analysis}}

\textbf{{What was being measured.}}
A final closing analysis was performed for PTCH1 and OLIG2 to determine whether their TFBS annotations could identify a specific transcription factor that might plausibly serve as a master regulator of the broader developmental bridge around the SIM2-centered methylation signal.

\textbf{{Methods.}}
Probe-level evidence for PTCH1 and OLIG2 was extracted from the focused candidate full-universe evidence table. Available regulatory fields were inspected, including TFBS labels, regulatory feature groups, DNase annotations, open-chromatin annotations, and shortlist presence in the original top-200 directional export. The TFBS field was examined specifically to determine whether it encoded named transcription factors or only genomic intervals.

\textbf{{Results.}}
{latex_escape(tfbs_note)}
The overall master-TF inference status for this closing step was: {latex_escape(master_status)}.

\textbf{{Tables written.}}
\texttt{{results/gse240704/ptch1\_olig2\_tfbs\_closing\_analysis/ptch1\_olig2\_full\_details.csv}}\\
\texttt{{results/gse240704/ptch1\_olig2\_tfbs\_closing\_analysis/ptch1\_olig2\_combined\_summary.csv}}\\
\texttt{{results/gse240704/ptch1\_olig2\_tfbs\_closing\_analysis/tables\_tex/table\_ptch1\_olig2\_tfbs\_closing.tex}}\\
\texttt{{results/gse240704/ptch1\_olig2\_tfbs\_closing\_analysis/ptch1\_olig2\_closing\_summary.json}}

\textbf{{Figures.}}
\texttt{{results/gse240704/ptch1\_olig2\_tfbs\_closing\_analysis/plots/ptch1\_olig2\_best\_effect.png}}

\textbf{{Interpretation.}}
This closing step tests whether the current methylation-annotation layer is sufficient to nominate a specific transcription factor as the master of the bridge cluster. If TFBS fields contain only genomic intervals and not named factors, then a master-TF claim is not earned from this dataset branch alone. In that case, the correct conclusion is that the branch supports a dominant SIM2-centered local methylation signature and a weaker distributed developmental backdrop, but not a specific transcription-factor master regulator.

\textbf{{Next steps.}}
This dataset branch can now be closed at the interpretive level. Any stronger master-regulator inference would require external TF-motif annotation, richer regulatory resources, or paired expression and chromatin data rather than the current methylation annotation fields alone.
""".strip() + "\n"

    outpath.write_text(tex, encoding="utf-8")


def main() -> None:
    ensure_dir(OUTDIR)
    ensure_dir(PLOTDIR)
    ensure_dir(TABLETEXDIR)
    ensure_dir(LOGBOOKDIR)

    if not FOCUSED_HITS_TABLE.exists():
        raise FileNotFoundError(f"Missing focused hits table: {FOCUSED_HITS_TABLE}")
    if not FOCUSED_SUMMARY_TABLE.exists():
        raise FileNotFoundError(f"Missing focused summary table: {FOCUSED_SUMMARY_TABLE}")
    if not TOP200_INTERSECTION_TABLE.exists():
        raise FileNotFoundError(f"Missing top200 intersection table: {TOP200_INTERSECTION_TABLE}")

    hits_df = pd.read_csv(FOCUSED_HITS_TABLE)
    summary_df = pd.read_csv(FOCUSED_SUMMARY_TABLE)
    top200_df = pd.read_csv(TOP200_INTERSECTION_TABLE)

    ptch1_olig2_details = hits_df[
        hits_df["candidate_gene"].astype(str).str.upper().isin(TARGET_GENES)
    ].copy()

    details_csv = OUTDIR / "ptch1_olig2_full_details.csv"
    ptch1_olig2_details.to_csv(details_csv, index=False)

    combined_summary_df = build_combined_summary(ptch1_olig2_details, top200_df)
    combined_summary_csv = OUTDIR / "ptch1_olig2_combined_summary.csv"
    combined_summary_df.to_csv(combined_summary_csv, index=False)

    tfbs_ptch1 = summarize_tfbs_identity(ptch1_olig2_details, "PTCH1")
    tfbs_olig2 = summarize_tfbs_identity(ptch1_olig2_details, "OLIG2")

    tfbs_identity_rows = pd.DataFrame([tfbs_ptch1, tfbs_olig2])
    tfbs_identity_csv = OUTDIR / "ptch1_olig2_tfbs_identity_summary.csv"
    tfbs_identity_rows.to_csv(tfbs_identity_csv, index=False)

    named_candidates = []
    for d in [tfbs_ptch1, tfbs_olig2]:
        named_candidates.extend(d["named_tf_candidates"])
    named_candidates = sorted(set(named_candidates))

    if named_candidates:
        overall_status = "named_tf_candidates_present"
        tfbs_statement = (
            "Named TF candidates were present in the TFBS annotation fields, "
            "so a tentative master-regulator discussion could in principle proceed."
        )
    else:
        all_coord_flags = [
            d["all_tfbs_entries_are_coordinates"]
            for d in [tfbs_ptch1, tfbs_olig2]
            if d["n_tfbs_annotated_probes"] > 0
        ]
        if all_coord_flags and all(all_coord_flags):
            overall_status = "tfbs_intervals_without_named_tf_identity"
            tfbs_statement = (
                "The PTCH1 and OLIG2 TFBS fields contained interval-style coordinate labels rather than named transcription-factor identities, "
                "so this dataset branch does not support nomination of a specific master TF from the current annotation layer."
            )
        else:
            overall_status = "no_resolved_master_tf_identity"
            tfbs_statement = (
                "No resolved named transcription-factor identity could be extracted from the PTCH1 and OLIG2 TFBS annotation fields."
            )

    evidence_tex = TABLETEXDIR / "table_ptch1_olig2_tfbs_closing.tex"
    write_summary_tex(combined_summary_df, evidence_tex)

    plot_path = PLOTDIR / "ptch1_olig2_best_effect.png"
    write_effect_plot(combined_summary_df, plot_path)

    matched_background_note = None
    if MATCHED_BACKGROUND_TABLE.exists():
        mb = pd.read_csv(MATCHED_BACKGROUND_TABLE)
        if not mb.empty and "feature" in mb.columns and "empirical_p_enrichment" in mb.columns:
            mb = mb.sort_values("empirical_p_enrichment", ascending=True).reset_index(drop=True)
            top_feat = mb.iloc[0]["feature"]
            top_p = mb.iloc[0]["empirical_p_enrichment"]
            matched_background_note = f"Strongest matched-background feature was {top_feat} with empirical p={top_p:.4f}."

    summary_json = {
        "dataset": DATASET,
        "target_genes": TARGET_GENES,
        "master_tf_inference_overall": overall_status,
        "tfbs_identity_statement": tfbs_statement,
        "named_tf_candidates": named_candidates,
        "ptch1_tfbs_identity": tfbs_ptch1,
        "olig2_tfbs_identity": tfbs_olig2,
        "genes_present_in_top200": combined_summary_df.loc[
            combined_summary_df["present_in_top200"], "candidate_gene"
        ].tolist(),
        "matched_background_note": matched_background_note,
        "outputs": {
            "details_csv": str(details_csv),
            "combined_summary_csv": str(combined_summary_csv),
            "tfbs_identity_csv": str(tfbs_identity_csv),
            "summary_tex": str(evidence_tex),
            "plot": str(plot_path),
        },
    }

    summary_json_path = OUTDIR / "ptch1_olig2_closing_summary.json"
    with open(summary_json_path, "w", encoding="utf-8") as f:
        json.dump(summary_json, f, indent=2)

    logbook_path = LOGBOOKDIR / "entry_35_ptch1_olig2_tfbs_closing_analysis.tex"
    write_logbook_entry(summary_json, logbook_path)

    print(json.dumps(summary_json, indent=2))
    print(f"[ok] wrote logbook entry to {logbook_path}")


if __name__ == "__main__":
    main()
