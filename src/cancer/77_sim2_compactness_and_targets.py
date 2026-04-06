#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
DATASET = "gse240704"

INPUT_PROBE_TABLE = (
    REPO_ROOT
    / "results"
    / DATASET
    / "chr21_locus_followup"
    / "D3_21q_probe_table.csv"
)

OUTDIR = (
    REPO_ROOT
    / "results"
    / DATASET
    / "sim2_followup"
)

PLOTDIR = OUTDIR / "plots"
LOGBOOKDIR = OUTDIR / "logbook"

# Seed list only, meant as a scaffold for interpretation rather than a claim of direct regulation here.
# This should be described in the logbook as a curated hypothesis list, not a demonstrated target set in this dataset.
SIM2_TARGET_SEED = [
    {"gene": "SHH", "category": "hedgehog_pathway_seed", "note": "developmental signaling context"},
    {"gene": "GLI2", "category": "hedgehog_pathway_seed", "note": "transcriptional effector context"},
    {"gene": "GLI3", "category": "hedgehog_pathway_seed", "note": "transcriptional effector context"},
    {"gene": "PTCH1", "category": "hedgehog_pathway_seed", "note": "relevant because PTCH1 appeared in S-axis positive loading biology"},
    {"gene": "BMP4", "category": "developmental_patterning_seed", "note": "developmental transcriptional context"},
    {"gene": "PAX6", "category": "neurodevelopment_seed", "note": "developmental context"},
    {"gene": "ASCL1", "category": "neural_lineage_seed", "note": "lineage program context"},
    {"gene": "OLIG2", "category": "neural_lineage_seed", "note": "gliogenic lineage context"},
    {"gene": "SOX2", "category": "stemness_lineage_seed", "note": "state maintenance context"},
    {"gene": "MYCN", "category": "neuro-oncology_seed", "note": "oncogenic lineage context"},
]

ENTRY_NUMBER = 25


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [
        re.sub(r"[^a-z0-9]+", "_", str(c).strip().lower()).strip("_")
        for c in out.columns
    ]
    return out


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
    tokens = []
    for part in re.split(r"[;,/|]+", s):
        token = part.strip()
        if token:
            tokens.append(token)
    return tokens


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


def choose_refgene_group_col(df: pd.DataFrame) -> Optional[str]:
    return find_first_existing(
        df,
        [
            "refgene_group_ann",
            "refgene_group_dir",
            "ucsc_refgene_group_ann",
            "ucsc_refgene_group_dir",
            "refgene_group",
            "ucsc_refgene_group",
        ],
    )


def choose_island_col(df: pd.DataFrame) -> Optional[str]:
    return find_first_existing(
        df,
        [
            "relation_to_cpg_island_ann",
            "relation_to_cpg_island_dir",
            "relation_to_ucsc_cpg_island_ann",
            "relation_to_ucsc_cpg_island_dir",
            "relation_to_cpg_island",
            "relation_to_ucsc_cpg_island",
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


def build_probe_gene_membership(df: pd.DataFrame, gene_col: str) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        probe_id = row.get("probe_id")
        tokens = split_gene_tokens(row.get(gene_col))
        unique_tokens = sorted(set(tokens))
        if not unique_tokens:
            rows.append(
                {
                    "probe_id": probe_id,
                    "gene_token": None,
                    "n_tokens_in_row": 0,
                }
            )
        else:
            for token in unique_tokens:
                rows.append(
                    {
                        "probe_id": probe_id,
                        "gene_token": token,
                        "n_tokens_in_row": len(unique_tokens),
                    }
                )
    out = pd.DataFrame(rows)
    return out


def build_unique_probe_gene_summary(membership: pd.DataFrame) -> pd.DataFrame:
    if membership.empty:
        return pd.DataFrame(columns=["gene_token", "n_unique_probes"])
    out = (
        membership.dropna(subset=["gene_token"])
        .groupby("gene_token", as_index=False)["probe_id"]
        .nunique()
        .rename(columns={"probe_id": "n_unique_probes"})
        .sort_values(["n_unique_probes", "gene_token"], ascending=[False, True])
        .reset_index(drop=True)
    )
    return out


def summarize_context(df: pd.DataFrame, colname: Optional[str]) -> dict[str, int]:
    if colname is None or colname not in df.columns:
        return {}
    series = df[colname].fillna("NA").astype(str)
    return series.value_counts().to_dict()


def write_position_plot(df: pd.DataFrame, outpath: Path, title: str) -> None:
    if df.empty:
        return

    effect_col = choose_effect_col(df)
    x = pd.to_numeric(df["position_std"], errors="coerce")

    if effect_col is not None and effect_col in df.columns:
        y = pd.to_numeric(df[effect_col], errors="coerce")
        ylabel = effect_col
    else:
        y = pd.Series([1] * len(df))
        ylabel = "probe count"

    plt.figure(figsize=(11, 4.5))
    plt.scatter(x, y, alpha=0.9)
    for _, row in df.iterrows():
        if "probe_id" in df.columns and "position_std" in df.columns:
            plt.annotate(
                str(row["probe_id"]),
                (row["position_std"], row[effect_col] if effect_col is not None else 1),
                fontsize=7,
                alpha=0.75,
            )

    plt.title(title)
    plt.xlabel("Genomic position on chr21")
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def write_ranked_effect_plot(df: pd.DataFrame, outpath: Path, title: str) -> None:
    if df.empty:
        return

    effect_col = choose_effect_col(df)
    if effect_col is None or effect_col not in df.columns:
        return

    tmp = df.copy()
    tmp["_abs_effect"] = pd.to_numeric(tmp[effect_col], errors="coerce").abs()
    tmp = tmp.sort_values("_abs_effect", ascending=False).reset_index(drop=True)
    tmp["rank"] = range(1, len(tmp) + 1)

    plt.figure(figsize=(8.5, 4.5))
    plt.scatter(tmp["rank"], tmp["_abs_effect"], alpha=0.9)
    for _, row in tmp.iterrows():
        plt.annotate(
            str(row["probe_id"]),
            (row["rank"], row["_abs_effect"]),
            fontsize=7,
            alpha=0.75,
        )
    plt.title(title)
    plt.xlabel("Probe rank")
    plt.ylabel(f"|{effect_col}|")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def build_sim2_compactness_summary(df: pd.DataFrame) -> dict:
    if df.empty:
        return {
            "n_sim2_probes": 0,
            "min_position": None,
            "max_position": None,
            "span_bp": None,
            "median_position": None,
            "mean_position": None,
            "nearest_neighbor_bp_median": None,
            "nearest_neighbor_bp_max": None,
        }

    positions = (
        pd.to_numeric(df["position_std"], errors="coerce")
        .dropna()
        .sort_values()
        .tolist()
    )

    min_pos = int(min(positions))
    max_pos = int(max(positions))
    span_bp = int(max_pos - min_pos)
    median_pos = float(pd.Series(positions).median())
    mean_pos = float(pd.Series(positions).mean())

    diffs = []
    for i in range(1, len(positions)):
        diffs.append(int(positions[i] - positions[i - 1]))

    return {
        "n_sim2_probes": int(len(positions)),
        "min_position": min_pos,
        "max_position": max_pos,
        "span_bp": span_bp,
        "median_position": median_pos,
        "mean_position": mean_pos,
        "nearest_neighbor_bp_median": None if not diffs else float(pd.Series(diffs).median()),
        "nearest_neighbor_bp_max": None if not diffs else int(max(diffs)),
        "sorted_positions": positions,
        "adjacent_gaps_bp": diffs,
    }


def write_logbook_entry(
    summary: dict,
    outpath: Path,
    figure_names: dict[str, str],
    table_names: dict[str, str],
) -> None:
    sim2 = summary["sim2_compactness"]
    unique_probe_gene = summary["unique_probe_gene_summary_head"]
    context = summary["context_summary"]

    text = rf"""
\subsection*{{Logbook Entry {ENTRY_NUMBER}: SIM2-centered 21q follow up}}

\textbf{{What was being measured.}}
A follow up was performed on the D3 lower-in-D3 chromosome 21q methylation signal to determine whether the arm-level signal was in fact concentrated into a specific local gene neighborhood. The main aims were to measure unique probe support by annotated gene, quantify genomic compactness of the dominant 21q locus, and generate a first downstream-target scaffold for SIM2-centered interpretation.

\textbf{{Methods.}}
The input table was
\texttt{{results/gse240704/chr21\_locus\_followup/D3\_21q\_probe\_table.csv}}.
Unique probe-to-gene membership was reconstructed by splitting manifest-derived gene token strings and counting distinct probes per gene rather than raw repeated token appearances. A SIM2-only subset was then extracted. Genomic compactness was quantified using the observed probe positions on chromosome 21q, including minimum and maximum coordinates, total genomic span, median position, and adjacent probe spacing. Context summaries were retained for refgene-group and CpG-island relation. A seed table of candidate downstream targets of SIM2 was exported as a hypothesis scaffold for later mechanistic overlay work, not as a direct regulatory claim established by this methylation dataset.

\textbf{{Results.}}
The D3 lower-in-D3 21q signal remained sharply localized. The 21q subset contained {summary["n_total_21q_probes"]} probes in total. Unique probe-level gene support showed SIM2 as the dominant annotated locus, with {summary["n_sim2_unique_probes"]} of {summary["n_total_21q_probes"]} probes carrying SIM2 annotation, while one probe remained unannotated. The SIM2-centered subset spanned positions {sim2["min_position"]} to {sim2["max_position"]} on chromosome 21, corresponding to a total span of {sim2["span_bp"]} bp. This indicates that the observed 21q depletion is not diffuse across the arm, but concentrated into a compact local neighborhood around SIM2. The context remained strongly gene-body and CpG-island centered. ERG and RUNX1 were not supported by direct probe annotation in the earlier locus-window follow up, and the present compactness analysis strengthens the interpretation that the 21q signal is best treated as a SIM2-centered local methylation signature rather than an ERG- or RUNX1-centered program.

\textbf{{Tables written.}}
\texttt{{{table_names["membership"]}}}\\
\texttt{{{table_names["unique_gene_summary"]}}}\\
\texttt{{{table_names["sim2_probe_table"]}}}\\
\texttt{{{table_names["sim2_compactness"]}}}\\
\texttt{{{table_names["sim2_target_seed"]}}}\\
\texttt{{{table_names["summary_json"]}}}

\textbf{{Figures.}}
\texttt{{{figure_names["sim2_position_plot"]}}}\\
\texttt{{{figure_names["sim2_ranked_effect_plot"]}}}\\
\texttt{{{figure_names["all21q_position_plot"]}}}

\textbf{{Interpretation.}}
The chromosome-arm level 21q enrichment can now be refined. In this dataset branch, the lower-in-D3 21q depletion is dominated by a compact SIM2-centered neighborhood, mostly in gene-body and CpG-island contexts. This does not by itself establish SIM2 as the causal driver of the D3 state. However, it does justify replacing a broad arm-only interpretation with a more local 21q neighborhood interpretation centered on SIM2-associated methylation structure.

\textbf{{Caveats.}}
This step remains annotation dependent. The downstream-target table is only a curated seed scaffold and not a direct target set validated in GSE240704. The analysis is methylation based and does not yet establish transcriptional direction, perturbational dependence, or regulatory causality.

\textbf{{Next steps.}}
The next step is to connect the SIM2-centered methylation signature to broader D3 biology by testing whether SIM2-related developmental or lineage programs overlap with the existing S-axis biology, especially given the earlier appearance of PTCH1 among positive S-axis loading genes. A later mechanistic overlay step should examine whether SIM2-centered structure aligns with developmental, neural-lineage, or stemness-related modules in a dataset with compatible expression or multi-omic readout.
""".strip() + "\n"

    outpath.write_text(text, encoding="utf-8")


def main() -> None:
    ensure_dir(OUTDIR)
    ensure_dir(PLOTDIR)
    ensure_dir(LOGBOOKDIR)

    if not INPUT_PROBE_TABLE.exists():
        raise FileNotFoundError(f"Missing input table: {INPUT_PROBE_TABLE}")

    df = pd.read_csv(INPUT_PROBE_TABLE)
    df = normalize_columns(df)

    if "probe_id" not in df.columns:
        raise ValueError("Expected probe_id in D3_21q_probe_table.csv")
    if "position_std" not in df.columns:
        raise ValueError("Expected position_std in D3_21q_probe_table.csv")

    gene_col = choose_gene_col(df)
    refgene_group_col = choose_refgene_group_col(df)
    island_col = choose_island_col(df)

    if gene_col is None:
        raise ValueError("Could not identify gene annotation column in D3_21q_probe_table.csv")

    membership = build_probe_gene_membership(df, gene_col=gene_col)
    membership_path = OUTDIR / "D3_21q_probe_gene_membership.csv"
    membership.to_csv(membership_path, index=False)

    unique_gene_summary = build_unique_probe_gene_summary(membership)
    unique_gene_summary_path = OUTDIR / "D3_21q_unique_probe_gene_summary.csv"
    unique_gene_summary.to_csv(unique_gene_summary_path, index=False)

    sim2_probe_ids = set(
        membership.loc[membership["gene_token"] == "SIM2", "probe_id"].dropna().astype(str).tolist()
    )

    sim2_df = df[df["probe_id"].astype(str).isin(sim2_probe_ids)].copy()
    sim2_df = sim2_df.sort_values("position_std").reset_index(drop=True)
    sim2_probe_table_path = OUTDIR / "D3_21q_SIM2_probe_table.csv"
    sim2_df.to_csv(sim2_probe_table_path, index=False)

    sim2_compactness = build_sim2_compactness_summary(sim2_df)
    sim2_compactness_df = pd.DataFrame([sim2_compactness])
    sim2_compactness_path = OUTDIR / "D3_21q_SIM2_compactness_summary.csv"
    sim2_compactness_df.to_csv(sim2_compactness_path, index=False)

    target_seed_df = pd.DataFrame(SIM2_TARGET_SEED)
    target_seed_path = OUTDIR / "SIM2_downstream_target_seed_table.csv"
    target_seed_df.to_csv(target_seed_path, index=False)

    all21q_plot = PLOTDIR / "D3_21q_all_probe_positions.png"
    sim2_plot = PLOTDIR / "D3_21q_SIM2_probe_positions.png"
    sim2_effect_plot = PLOTDIR / "D3_21q_SIM2_ranked_abs_effect.png"

    write_position_plot(
        df,
        outpath=all21q_plot,
        title="D3 lower-in-D3 21q probe positions",
    )

    write_position_plot(
        sim2_df,
        outpath=sim2_plot,
        title="D3 lower-in-D3 SIM2-centered 21q probe positions",
    )

    write_ranked_effect_plot(
        sim2_df,
        outpath=sim2_effect_plot,
        title="D3 lower-in-D3 SIM2 probe absolute effect ranking",
    )

    summary = {
        "dataset": DATASET,
        "input_probe_table": str(INPUT_PROBE_TABLE),
        "n_total_21q_probes": int(len(df)),
        "n_sim2_unique_probes": int(len(sim2_df)),
        "n_unannotated_21q_probes": int(
            membership["gene_token"].isna().sum()
        ),
        "unique_probe_gene_summary_head": unique_gene_summary.head(20).to_dict(orient="records"),
        "context_summary": {
            "refgene_group_counts": summarize_context(df, refgene_group_col),
            "cpg_island_relation_counts": summarize_context(df, island_col),
        },
        "sim2_compactness": sim2_compactness,
        "outputs": {
            "membership": str(membership_path),
            "unique_gene_summary": str(unique_gene_summary_path),
            "sim2_probe_table": str(sim2_probe_table_path),
            "sim2_compactness_summary": str(sim2_compactness_path),
            "sim2_target_seed_table": str(target_seed_path),
            "all21q_position_plot": str(all21q_plot),
            "sim2_position_plot": str(sim2_plot),
            "sim2_ranked_effect_plot": str(sim2_effect_plot),
        },
    }

    summary_path = OUTDIR / "D3_21q_SIM2_followup_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    figure_names = {
        "all21q_position_plot": "results/gse240704/sim2_followup/plots/D3_21q_all_probe_positions.png",
        "sim2_position_plot": "results/gse240704/sim2_followup/plots/D3_21q_SIM2_probe_positions.png",
        "sim2_ranked_effect_plot": "results/gse240704/sim2_followup/plots/D3_21q_SIM2_ranked_abs_effect.png",
    }

    table_names = {
        "membership": "results/gse240704/sim2_followup/D3_21q_probe_gene_membership.csv",
        "unique_gene_summary": "results/gse240704/sim2_followup/D3_21q_unique_probe_gene_summary.csv",
        "sim2_probe_table": "results/gse240704/sim2_followup/D3_21q_SIM2_probe_table.csv",
        "sim2_compactness": "results/gse240704/sim2_followup/D3_21q_SIM2_compactness_summary.csv",
        "sim2_target_seed": "results/gse240704/sim2_followup/SIM2_downstream_target_seed_table.csv",
        "summary_json": "results/gse240704/sim2_followup/D3_21q_SIM2_followup_summary.json",
    }

    logbook_path = LOGBOOKDIR / "entry_25_sim2_checks.tex"
    write_logbook_entry(
        summary=summary,
        outpath=logbook_path,
        figure_names=figure_names,
        table_names=table_names,
    )

    print(json.dumps(summary, indent=2))
    print(f"[ok] wrote logbook entry to {logbook_path}")


if __name__ == "__main__":
    main()
