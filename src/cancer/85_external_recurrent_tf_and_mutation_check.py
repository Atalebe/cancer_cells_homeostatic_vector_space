#!/usr/bin/env python3

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Optional

import matplotlib.pyplot as plt
import pandas as pd


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
DATASET = "gse240704"

EXTERNAL_SUMMARY_JSON = (
    REPO_ROOT
    / "results"
    / DATASET
    / "external_tfbs_overlap"
    / "summary"
    / "external_tfbs_overlap_summary.json"
)

INTERVAL_LABELS = [
    "PTCH1_tfbs_interval",
    "OLIG2_tfbs_interval",
    "SIM2_core_cluster",
    "PAX6_best_probe_window",
]

OUTDIR = (
    REPO_ROOT
    / "results"
    / DATASET
    / "external_tf_recurrence_and_mutation_check"
)

PLOTDIR = OUTDIR / "plots"
TABLETEXDIR = OUTDIR / "tables_tex"
LOGBOOKDIR = OUTDIR / "logbook"

ENTRY_NUMBER = 36


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


def load_external_summary() -> list[dict[str, Any]]:
    if not EXTERNAL_SUMMARY_JSON.exists():
        raise FileNotFoundError(f"Missing external summary json: {EXTERNAL_SUMMARY_JSON}")
    with open(EXTERNAL_SUMMARY_JSON, "r", encoding="utf-8") as f:
        obj = json.load(f)
    if not isinstance(obj, list):
        raise ValueError("Expected external summary JSON to be a list.")
    return obj


def build_recurrent_tf_table(records: list[dict[str, Any]]) -> pd.DataFrame:
    rows = []

    for rec in records:
        if rec.get("source") != "ucsc":
            continue
        label = rec.get("label")
        tf_map = rec.get("top_tf_like_values", {}) or {}
        if label not in INTERVAL_LABELS:
            continue
        for tf_name, count in tf_map.items():
            rows.append(
                {
                    "interval_label": label,
                    "tf_name": str(tf_name),
                    "count_in_interval": int(count),
                }
            )

    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(
            columns=[
                "tf_name",
                "n_intervals_present",
                "total_count_across_intervals",
                "intervals_present",
                "ptch1_count",
                "olig2_count",
                "sim2_count",
                "pax6_count",
            ]
        )

    summary_rows = []
    for tf_name, sub in df.groupby("tf_name", sort=True):
        interval_counts = dict(zip(sub["interval_label"], sub["count_in_interval"]))
        intervals_present = sorted(sub["interval_label"].tolist())
        summary_rows.append(
            {
                "tf_name": tf_name,
                "n_intervals_present": int(len(intervals_present)),
                "total_count_across_intervals": int(sub["count_in_interval"].sum()),
                "intervals_present": "; ".join(intervals_present),
                "ptch1_count": int(interval_counts.get("PTCH1_tfbs_interval", 0)),
                "olig2_count": int(interval_counts.get("OLIG2_tfbs_interval", 0)),
                "sim2_count": int(interval_counts.get("SIM2_core_cluster", 0)),
                "pax6_count": int(interval_counts.get("PAX6_best_probe_window", 0)),
            }
        )

    out = pd.DataFrame(summary_rows)
    out = out.sort_values(
        ["n_intervals_present", "total_count_across_intervals", "tf_name"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    return out


def find_candidate_mutation_files() -> list[Path]:
    search_roots = [
        REPO_ROOT / "data",
        REPO_ROOT / "results" / DATASET,
    ]

    include_name_patterns = [
        "mut",
        "mutation",
        "variant",
        "maf",
        "vcf",
        "snv",
        "indel",
        "oncoprint",
        "oncoplot",
    ]
    exts = {".csv", ".tsv", ".txt", ".maf", ".vcf"}

    found = []
    for root in search_roots:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            if path.suffix.lower() not in exts:
                continue
            name_l = path.name.lower()
            if any(pat in name_l for pat in include_name_patterns):
                found.append(path)

    return sorted(set(found))


def guess_delimiter(path: Path) -> str:
    if path.suffix.lower() in {".tsv", ".maf", ".vcf"}:
        return "\t"
    return ","


def read_tabular(path: Path) -> Optional[pd.DataFrame]:
    try:
        return pd.read_csv(path, sep=guess_delimiter(path), low_memory=False, comment="#")
    except Exception:
        return None


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [
        re.sub(r"[^a-z0-9]+", "_", str(c).strip().lower()).strip("_")
        for c in out.columns
    ]
    return out


def find_gene_column(df: pd.DataFrame) -> Optional[str]:
    candidates = [
        "hugo_symbol",
        "gene",
        "gene_symbol",
        "symbol",
        "gene_name",
        "tumor_gene",
        "genes",
    ]
    for c in candidates:
        if c in df.columns:
            return c
    return None


def build_mutation_hits(recurrent_tfs: list[str], mutation_files: list[Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    file_rows = []
    hit_rows = []

    recurrent_set = {x.upper() for x in recurrent_tfs}

    for path in mutation_files:
        df = read_tabular(path)
        if df is None or df.empty:
            file_rows.append(
                {
                    "path": str(path),
                    "status": "failed_or_empty",
                    "n_rows": None,
                    "gene_column": None,
                    "n_hits": 0,
                }
            )
            continue

        df = normalize_columns(df)
        gene_col = find_gene_column(df)

        if gene_col is None:
            file_rows.append(
                {
                    "path": str(path),
                    "status": "no_gene_column",
                    "n_rows": int(len(df)),
                    "gene_column": None,
                    "n_hits": 0,
                }
            )
            continue

        gene_series = df[gene_col].fillna("").astype(str).str.upper()
        mask = gene_series.isin(recurrent_set)
        sub = df.loc[mask].copy()

        file_rows.append(
            {
                "path": str(path),
                "status": "scanned",
                "n_rows": int(len(df)),
                "gene_column": gene_col,
                "n_hits": int(len(sub)),
            }
        )

        if not sub.empty:
            sub["matched_tf"] = gene_series.loc[sub.index].values
            sub["source_file"] = str(path)
            hit_rows.append(sub)

    file_df = pd.DataFrame(file_rows).sort_values(["n_hits", "path"], ascending=[False, True]).reset_index(drop=True)
    if hit_rows:
        hits_df = pd.concat(hit_rows, axis=0, ignore_index=True)
    else:
        hits_df = pd.DataFrame()

    return file_df, hits_df


def build_mutation_summary(hits_df: pd.DataFrame, recurrent_tfs: list[str]) -> pd.DataFrame:
    rows = []
    recurrent_set = [x.upper() for x in recurrent_tfs]

    if hits_df.empty:
        return pd.DataFrame(
            {
                "tf_name": recurrent_set,
                "n_rows_detected": [0] * len(recurrent_set),
            }
        )

    for tf in recurrent_set:
        sub = hits_df[hits_df["matched_tf"] == tf]
        rows.append(
            {
                "tf_name": tf,
                "n_rows_detected": int(len(sub)),
                "source_files": "; ".join(sorted(set(sub["source_file"].astype(str).tolist()))),
            }
        )

    out = pd.DataFrame(rows).sort_values(["n_rows_detected", "tf_name"], ascending=[False, True]).reset_index(drop=True)
    return out


def write_recurrent_tf_tex(df: pd.DataFrame, outpath: Path) -> None:
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Recurrent external TF names across PTCH1, OLIG2, SIM2, and PAX6 interval overlaps from the UCSC hg19 ENCODE TFBS clustered track.}")
    lines.append(r"\label{tab:external_recurrent_tf}")
    lines.append(r"\begin{tabular}{lccc}")
    lines.append(r"\hline")
    lines.append(r"TF & Intervals & Total count & Intervals present \\")
    lines.append(r"\hline")

    for _, row in df.iterrows():
        tf = latex_escape(row["tf_name"])
        n_intervals = latex_escape(row["n_intervals_present"])
        total = latex_escape(row["total_count_across_intervals"])
        intervals = latex_escape(row["intervals_present"])
        lines.append(f"{tf} & {n_intervals} & {total} & {intervals} \\\\")

    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    outpath.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_mutation_summary_tex(df: pd.DataFrame, outpath: Path) -> None:
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Repo-local mutation cross-reference for recurrent external TF candidates. Presence depends on whether any mutation or variant table exists in the current dataset branch or repo structure.}")
    lines.append(r"\label{tab:external_recurrent_tf_mutation_check}")
    lines.append(r"\begin{tabular}{lc}")
    lines.append(r"\hline")
    lines.append(r"TF & Mutation-table rows detected \\")
    lines.append(r"\hline")

    for _, row in df.iterrows():
        tf = latex_escape(row["tf_name"])
        n = latex_escape(row["n_rows_detected"])
        lines.append(f"{tf} & {n} \\\\")

    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    outpath.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_plot(df: pd.DataFrame, outpath: Path) -> None:
    if df.empty:
        return
    tmp = df.copy().head(12)
    tmp["rank"] = range(1, len(tmp) + 1)

    plt.figure(figsize=(10, 4.8))
    plt.scatter(tmp["rank"], tmp["n_intervals_present"], alpha=0.9)
    for _, row in tmp.iterrows():
        plt.annotate(
            str(row["tf_name"]),
            (row["rank"], row["n_intervals_present"]),
            fontsize=8,
            alpha=0.8,
        )
    plt.title("Recurrent external TF names across queried regulatory intervals")
    plt.xlabel("TF rank")
    plt.ylabel("Number of queried intervals present")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def write_logbook_entry(summary_json: dict, outpath: Path) -> None:
    top_tfs = ", ".join(summary_json["top_recurrent_tfs"]) if summary_json["top_recurrent_tfs"] else "none"
    mutation_status = summary_json["mutation_check_status"]

    tex = rf"""
\subsection*{{Logbook Entry {ENTRY_NUMBER}: External recurrent TF closing pass and mutation cross-reference}}

\textbf{{What was being measured.}}
A final external closing pass was performed to identify recurrent transcription-factor names across the PTCH1, OLIG2, SIM2, and PAX6 regulatory intervals, and to determine whether those recurrent factors themselves could be cross-referenced against any mutation or variant tables present in the repo.

\textbf{{Methods.}}
The UCSC hg19 ENCODE TFBS clustered overlap summary was parsed across the four queried intervals. TF names were aggregated by interval, total count, and number of queried intervals in which they appeared. The recurrent TF candidates were then cross-referenced against any repo-local mutation-like tables discovered under the project data and results paths, including files with names suggestive of mutation, variant, MAF, VCF, SNV, or indel content. A mutation hit was defined as a row in such a table whose gene column matched one of the recurrent TF names.

\textbf{{Results.}}
The recurrent external TF candidates across the queried intervals were led by: {latex_escape(top_tfs)}.
Mutation cross-reference status for this dataset branch was: {latex_escape(mutation_status)}.

\textbf{{Interpretation.}}
This final pass does not by itself establish a master regulator, but it provides a cleaner external recurrence-based view of which regulatory factors repeatedly overlap the SIM2-centered and bridge-associated loci. Mutation cross-reference can strengthen mechanistic interest only if mutation data actually exist in the repo and show corresponding alterations in those same factors. Absence of mutation hits, or complete absence of mutation tables, should not be overinterpreted.

\textbf{{Tables written.}}
\texttt{{results/gse240704/external\_tf\_recurrence\_and\_mutation\_check/external\_recurrent\_tf\_summary.csv}}\\
\texttt{{results/gse240704/external\_tf\_recurrence\_and\_mutation\_check/mutation\_candidate\_files.csv}}\\
\texttt{{results/gse240704/external\_tf\_recurrence\_and\_mutation\_check/recurrent\_tf\_mutation\_summary.csv}}\\
\texttt{{results/gse240704/external\_tf\_recurrence\_and\_mutation\_check/tables\_tex/table\_external\_recurrent\_tf.tex}}\\
\texttt{{results/gse240704/external\_tf\_recurrence\_and\_mutation\_check/tables\_tex/table\_external\_recurrent\_tf\_mutation\_check.tex}}\\
\texttt{{results/gse240704/external\_tf\_recurrence\_and\_mutation\_check/external\_recurrent\_tf\_summary.json}}

\textbf{{Figures.}}
\texttt{{results/gse240704/external\_tf\_recurrence\_and\_mutation\_check/plots/external\_recurrent\_tf\_rank.png}}

\textbf{{Next steps.}}
Use the external recurrence results to finalize the mechanistic language for the dataset. If no mutation layer is present or no recurrent-factor mutations are detected, then the correct closing frame remains locus-centric and regulatory-context aware, but not mutation-anchored.
""".strip() + "\n"

    outpath.write_text(tex, encoding="utf-8")


def main() -> None:
    ensure_dir(OUTDIR)
    ensure_dir(PLOTDIR)
    ensure_dir(TABLETEXDIR)
    ensure_dir(LOGBOOKDIR)

    records = load_external_summary()
    recurrent_df = build_recurrent_tf_table(records)
    recurrent_csv = OUTDIR / "external_recurrent_tf_summary.csv"
    recurrent_df.to_csv(recurrent_csv, index=False)

    top_recurrent_tfs = recurrent_df.head(10)["tf_name"].tolist()

    mutation_files = find_candidate_mutation_files()
    mutation_files_df, mutation_hits_df = build_mutation_hits(top_recurrent_tfs, mutation_files)

    mutation_files_csv = OUTDIR / "mutation_candidate_files.csv"
    mutation_files_df.to_csv(mutation_files_csv, index=False)

    mutation_hits_csv = OUTDIR / "recurrent_tf_mutation_hits.csv"
    mutation_hits_df.to_csv(mutation_hits_csv, index=False)

    mutation_summary_df = build_mutation_summary(mutation_hits_df, top_recurrent_tfs)
    mutation_summary_csv = OUTDIR / "recurrent_tf_mutation_summary.csv"
    mutation_summary_df.to_csv(mutation_summary_csv, index=False)

    recurrent_tex = TABLETEXDIR / "table_external_recurrent_tf.tex"
    write_recurrent_tf_tex(recurrent_df.head(15), recurrent_tex)

    mutation_tex = TABLETEXDIR / "table_external_recurrent_tf_mutation_check.tex"
    write_mutation_summary_tex(mutation_summary_df, mutation_tex)

    plot_path = PLOTDIR / "external_recurrent_tf_rank.png"
    write_plot(recurrent_df, plot_path)

    if len(mutation_files_df) == 0:
        mutation_status = "no_mutation_like_tables_found_in_repo_paths"
    elif mutation_hits_df.empty:
        mutation_status = "mutation_like_tables_scanned_but_no_recurrent_tf_hits_found"
    else:
        mutation_status = "mutation_hits_found_for_one_or_more_recurrent_tfs"

    summary_json = {
        "dataset": DATASET,
        "top_recurrent_tfs": top_recurrent_tfs,
        "n_mutation_like_files_found": int(len(mutation_files_df)),
        "mutation_check_status": mutation_status,
        "outputs": {
            "recurrent_tf_csv": str(recurrent_csv),
            "mutation_candidate_files_csv": str(mutation_files_csv),
            "mutation_hits_csv": str(mutation_hits_csv),
            "mutation_summary_csv": str(mutation_summary_csv),
            "recurrent_tf_tex": str(recurrent_tex),
            "mutation_tex": str(mutation_tex),
            "plot": str(plot_path),
        },
    }

    summary_json_path = OUTDIR / "external_recurrent_tf_summary.json"
    with open(summary_json_path, "w", encoding="utf-8") as f:
        json.dump(summary_json, f, indent=2)

    logbook_path = LOGBOOKDIR / "entry_36_external_recurrent_tf_and_mutation_check.tex"
    write_logbook_entry(summary_json, logbook_path)

    print(json.dumps(summary_json, indent=2))
    print(f"[ok] wrote logbook entry to {logbook_path}")


def find_candidate_mutation_files() -> list[Path]:
    search_roots = [
        REPO_ROOT / "data",
        REPO_ROOT / "results" / DATASET,
    ]

    include_name_patterns = [
        "mut",
        "mutation",
        "variant",
        "maf",
        "vcf",
        "snv",
        "indel",
        "oncoprint",
        "oncoplot",
    ]
    exts = {".csv", ".tsv", ".txt", ".maf", ".vcf"}

    found = []
    for root in search_roots:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            if path.suffix.lower() not in exts:
                continue
            name_l = path.name.lower()
            if any(pat in name_l for pat in include_name_patterns):
                found.append(path)

    return sorted(set(found))


def guess_delimiter(path: Path) -> str:
    if path.suffix.lower() in {".tsv", ".maf", ".vcf"}:
        return "\t"
    return ","


def read_tabular(path: Path) -> Optional[pd.DataFrame]:
    try:
        return pd.read_csv(path, sep=guess_delimiter(path), low_memory=False, comment="#")
    except Exception:
        return None


def normalize_columns_generic(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [
        re.sub(r"[^a-z0-9]+", "_", str(c).strip().lower()).strip("_")
        for c in out.columns
    ]
    return out


def find_gene_column(df: pd.DataFrame) -> Optional[str]:
    candidates = [
        "hugo_symbol",
        "gene",
        "gene_symbol",
        "symbol",
        "gene_name",
        "tumor_gene",
        "genes",
    ]
    for c in candidates:
        if c in df.columns:
            return c
    return None


def build_mutation_hits(recurrent_tfs: list[str], mutation_files: list[Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    file_rows = []
    hit_rows = []
    recurrent_set = {x.upper() for x in recurrent_tfs}

    for path in mutation_files:
        df = read_tabular(path)
        if df is None or df.empty:
            file_rows.append(
                {
                    "path": str(path),
                    "status": "failed_or_empty",
                    "n_rows": None,
                    "gene_column": None,
                    "n_hits": 0,
                }
            )
            continue

        df = normalize_columns_generic(df)
        gene_col = find_gene_column(df)

        if gene_col is None:
            file_rows.append(
                {
                    "path": str(path),
                    "status": "no_gene_column",
                    "n_rows": int(len(df)),
                    "gene_column": None,
                    "n_hits": 0,
                }
            )
            continue

        gene_series = df[gene_col].fillna("").astype(str).str.upper()
        mask = gene_series.isin(recurrent_set)
        sub = df.loc[mask].copy()

        file_rows.append(
            {
                "path": str(path),
                "status": "scanned",
                "n_rows": int(len(df)),
                "gene_column": gene_col,
                "n_hits": int(len(sub)),
            }
        )

        if not sub.empty:
            sub["matched_tf"] = gene_series.loc[sub.index].values
            sub["source_file"] = str(path)
            hit_rows.append(sub)

    file_df = pd.DataFrame(file_rows)
    if file_df.empty:
        file_df = pd.DataFrame(columns=["path", "status", "n_rows", "gene_column", "n_hits"])
    else:
        file_df = file_df.sort_values(["n_hits", "path"], ascending=[False, True]).reset_index(drop=True)

    if hit_rows:
        hits_df = pd.concat(hit_rows, axis=0, ignore_index=True)
    else:
        hits_df = pd.DataFrame()

    return file_df, hits_df


def build_mutation_summary(hits_df: pd.DataFrame, recurrent_tfs: list[str]) -> pd.DataFrame:
    rows = []
    recurrent_set = [x.upper() for x in recurrent_tfs]

    if hits_df.empty:
        return pd.DataFrame(
            {
                "tf_name": recurrent_set,
                "n_rows_detected": [0] * len(recurrent_set),
            }
        )

    for tf in recurrent_set:
        sub = hits_df[hits_df["matched_tf"] == tf]
        rows.append(
            {
                "tf_name": tf,
                "n_rows_detected": int(len(sub)),
                "source_files": "; ".join(sorted(set(sub["source_file"].astype(str).tolist()))),
            }
        )

    out = pd.DataFrame(rows).sort_values(["n_rows_detected", "tf_name"], ascending=[False, True]).reset_index(drop=True)
    return out


if __name__ == "__main__":
    main()
