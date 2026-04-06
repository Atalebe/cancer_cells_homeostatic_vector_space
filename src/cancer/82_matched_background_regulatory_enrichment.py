#!/usr/bin/env python3

from __future__ import annotations

import json
import random
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

MANIFEST_TABLE = (
    REPO_ROOT
    / "results"
    / DATASET
    / "manual_manifest_annotation"
    / "gpl23976_annotation_from_manifest.parquet"
)

FOCUSED_HITS_TABLE = (
    REPO_ROOT
    / "results"
    / DATASET
    / "focused_candidate_evidence"
    / "focused_candidate_full_universe_hits.csv"
)

OUTDIR = (
    REPO_ROOT
    / "results"
    / DATASET
    / "matched_background_regulatory_enrichment"
)

PLOTDIR = OUTDIR / "plots"
TABLETEXDIR = OUTDIR / "tables_tex"
LOGBOOKDIR = OUTDIR / "logbook"

ENTRY_NUMBER = 33
SEED = 123
N_MATCHED_SETS = 2000

FOCUSED_GENES = ["PAX6", "GLI2", "BMP4", "PTCH1", "SHH", "OLIG2"]
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


def load_universe_with_manifest() -> pd.DataFrame:
    if not FULL_UNIVERSE_TABLE.exists():
        raise FileNotFoundError(f"Missing full universe table: {FULL_UNIVERSE_TABLE}")
    if not MANIFEST_TABLE.exists():
        raise FileNotFoundError(f"Missing manifest table: {MANIFEST_TABLE}")

    base = normalize_columns(pd.read_csv(FULL_UNIVERSE_TABLE))
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


def choose_cols(df: pd.DataFrame) -> dict[str, Optional[str]]:
    return {
        "probe_id": find_first_existing(df, ["probe_id", "ilmnid", "ilmn_id", "name", "id_ref", "id"]),
        "gene": find_first_existing(
            df,
            [
                "ucsc_refgene_name_ann",
                "ucsc_refgene_name_dir",
                "ucsc_refgene_name",
                "gene_symbol",
                "gene_symbols",
                "genes",
            ],
        ),
        "abs_delta": find_first_existing(df, ["abs_delta_median"]),
        "delta": find_first_existing(df, ["delta_median_a_minus_b"]),
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


def add_matching_features(df: pd.DataFrame) -> pd.DataFrame:
    cols = choose_cols(df)
    out = df.copy()

    if cols["abs_delta"] and cols["abs_delta"] in out.columns:
        out["abs_effect_for_match"] = pd.to_numeric(out[cols["abs_delta"]], errors="coerce")
    elif cols["delta"] and cols["delta"] in out.columns:
        out["abs_effect_for_match"] = pd.to_numeric(out[cols["delta"]], errors="coerce").abs()
    else:
        raise ValueError("Could not identify an effect-size column for matching.")

    if cols["refgene_group"] and cols["refgene_group"] in out.columns:
        out["refgene_group_match"] = out[cols["refgene_group"]].fillna("NA").astype(str)
    else:
        out["refgene_group_match"] = "NA"

    if cols["cpg_relation"] and cols["cpg_relation"] in out.columns:
        out["cpg_relation_match"] = out[cols["cpg_relation"]].fillna("NA").astype(str)
    else:
        out["cpg_relation_match"] = "NA"

    def binary_from_nonempty(colname: Optional[str]) -> pd.Series:
        if colname is None or colname not in out.columns:
            return pd.Series([False] * len(out), index=out.index)
        return out[colname].fillna("").astype(str).str.strip().ne("")

    out["has_dnase"] = binary_from_nonempty(cols["dnase"])
    out["has_tfbs"] = binary_from_nonempty(cols["tfbs"])
    out["has_openchromatin"] = binary_from_nonempty(cols["openchromatin"])
    out["has_phantom4"] = binary_from_nonempty(cols["phantom4"])

    if cols["reg_feature_name"] and cols["reg_feature_name"] in out.columns:
        reg_name = out[cols["reg_feature_name"]].fillna("").astype(str)
        out["has_enhancer_like"] = reg_name.str.contains("enhancer", case=False, na=False)
    else:
        out["has_enhancer_like"] = False

    if cols["reg_feature_group"] and cols["reg_feature_group"] in out.columns:
        reg_group = out[cols["reg_feature_group"]].fillna("").astype(str)
        out["has_unclassified_cell_type_specific"] = reg_group.eq("Unclassified_Cell_type_specific")
    else:
        out["has_unclassified_cell_type_specific"] = False

    return out


def load_focused_probe_ids() -> list[str]:
    if not FOCUSED_HITS_TABLE.exists():
        raise FileNotFoundError(f"Missing focused hits table: {FOCUSED_HITS_TABLE}")
    hits = pd.read_csv(FOCUSED_HITS_TABLE)
    if "probe_id" not in hits.columns:
        raise ValueError("focused hits table missing probe_id")
    return sorted(hits["probe_id"].dropna().astype(str).unique().tolist())


def make_gene_token_set(value: object) -> set[str]:
    return {tok.upper() for tok in split_gene_tokens(value)}


def build_background_pool(df: pd.DataFrame, focused_probe_ids: set[str]) -> pd.DataFrame:
    cols = choose_cols(df)
    if cols["probe_id"] is None:
        raise ValueError("Universe table missing probe_id column.")
    if cols["gene"] is None:
        raise ValueError("Universe table missing gene annotation column.")

    out = df.copy()
    out["probe_id_str"] = out[cols["probe_id"]].astype(str)

    gene_series = out[cols["gene"]].fillna("").astype(str)

    excluded_genes = {g.upper() for g in FOCUSED_GENES + [LOCAL_ANCHOR]}

    def overlaps_excluded_gene(s: str) -> bool:
        tokens = make_gene_token_set(s)
        return any(g in tokens for g in excluded_genes)

    out = out[~out["probe_id_str"].isin(focused_probe_ids)].copy()
    out = out[~gene_series.map(overlaps_excluded_gene)].copy()
    out = out.dropna(subset=["abs_effect_for_match"]).copy()
    return out.reset_index(drop=True)


def sample_one_matched_set(
    focused_df: pd.DataFrame,
    background_df: pd.DataFrame,
    rng: random.Random,
    effect_tol_start: float = 0.01,
    effect_tol_max: float = 0.08,
) -> pd.DataFrame:
    chosen_indices: set[int] = set()
    sampled_rows = []

    for _, row in focused_df.iterrows():
        tol = effect_tol_start
        sub = pd.DataFrame()

        while tol <= effect_tol_max:
            sub = background_df[
                (background_df["refgene_group_match"] == row["refgene_group_match"])
                & (background_df["cpg_relation_match"] == row["cpg_relation_match"])
                & (background_df["abs_effect_for_match"] >= row["abs_effect_for_match"] - tol)
                & (background_df["abs_effect_for_match"] <= row["abs_effect_for_match"] + tol)
            ].copy()
            if len(sub) > 0:
                break
            tol *= 2.0

        if len(sub) == 0:
            tol = effect_tol_start
            while tol <= effect_tol_max:
                sub = background_df[
                    (background_df["refgene_group_match"] == row["refgene_group_match"])
                    & (background_df["abs_effect_for_match"] >= row["abs_effect_for_match"] - tol)
                    & (background_df["abs_effect_for_match"] <= row["abs_effect_for_match"] + tol)
                ].copy()
                if len(sub) > 0:
                    break
                tol *= 2.0

        if len(sub) == 0:
            tol = effect_tol_start
            while tol <= effect_tol_max:
                sub = background_df[
                    (background_df["abs_effect_for_match"] >= row["abs_effect_for_match"] - tol)
                    & (background_df["abs_effect_for_match"] <= row["abs_effect_for_match"] + tol)
                ].copy()
                if len(sub) > 0:
                    break
                tol *= 2.0

        if len(sub) == 0:
            raise RuntimeError("Could not find matched background probe for one focused probe.")

        available_idx = [i for i in sub.index.tolist() if i not in chosen_indices]
        if not available_idx:
            available_idx = sub.index.tolist()

        pick_idx = rng.choice(available_idx)
        chosen_indices.add(pick_idx)
        sampled_rows.append(background_df.loc[pick_idx])

    return pd.DataFrame(sampled_rows).reset_index(drop=True)


def summarize_feature_counts(df: pd.DataFrame) -> dict[str, int]:
    features = [
        "has_dnase",
        "has_tfbs",
        "has_openchromatin",
        "has_phantom4",
        "has_enhancer_like",
        "has_unclassified_cell_type_specific",
    ]
    return {feat: int(df[feat].sum()) for feat in features}


def empirical_pvalue(obs: int, null_values: list[int]) -> float:
    if not null_values:
        return float("nan")
    ge = sum(1 for x in null_values if x >= obs)
    return (ge + 1.0) / (len(null_values) + 1.0)


def build_enrichment_summary(observed: dict[str, int], null_summaries: list[dict[str, int]], n_probes: int) -> pd.DataFrame:
    rows = []
    features = list(observed.keys())

    for feat in features:
        null_vals = [d[feat] for d in null_summaries]
        null_mean = sum(null_vals) / len(null_vals)
        null_median = float(pd.Series(null_vals).median())
        p_emp = empirical_pvalue(observed[feat], null_vals)

        obs_rate = observed[feat] / n_probes if n_probes else float("nan")
        null_rate_mean = null_mean / n_probes if n_probes else float("nan")

        fold = None
        if null_rate_mean > 0:
            fold = obs_rate / null_rate_mean

        rows.append(
            {
                "feature": feat,
                "observed_count": observed[feat],
                "observed_rate": obs_rate,
                "null_mean_count": null_mean,
                "null_median_count": null_median,
                "null_mean_rate": null_rate_mean,
                "empirical_p_enrichment": p_emp,
                "fold_vs_null_mean_rate": fold,
            }
        )

    out = pd.DataFrame(rows).sort_values(
        ["empirical_p_enrichment", "feature"],
        ascending=[True, True],
    ).reset_index(drop=True)
    return out


def write_summary_tex(df: pd.DataFrame, outpath: Path) -> None:
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Matched-background enrichment test for regulatory annotations among the focused developmental bridge probes. Background probes were matched to the focused set on methylation-effect neighborhood and probe-context features.}")
    lines.append(r"\label{tab:matched_background_regulatory_enrichment}")
    lines.append(r"\begin{tabular}{lcccc}")
    lines.append(r"\hline")
    lines.append(r"Feature & Observed & Null mean & Fold & Empirical $p$ \\")
    lines.append(r"\hline")

    for _, row in df.iterrows():
        feature = latex_escape(row["feature"])
        obs = latex_escape(row["observed_count"])
        null_mean = f"{float(row['null_mean_count']):.2f}"
        fold = "" if pd.isna(row["fold_vs_null_mean_rate"]) else f"{float(row['fold_vs_null_mean_rate']):.2f}"
        pval = "" if pd.isna(row["empirical_p_enrichment"]) else f"{float(row['empirical_p_enrichment']):.4f}"
        lines.append(f"{feature} & {obs} & {null_mean} & {fold} & {pval} \\\\")

    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    outpath.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_feature_plot(df: pd.DataFrame, outpath: Path) -> None:
    tmp = df.copy()
    tmp = tmp.sort_values("empirical_p_enrichment", ascending=True).reset_index(drop=True)
    tmp["rank"] = range(1, len(tmp) + 1)

    plt.figure(figsize=(10, 4.8))
    plt.scatter(tmp["rank"], tmp["fold_vs_null_mean_rate"], alpha=0.9)
    for _, row in tmp.iterrows():
        y = row["fold_vs_null_mean_rate"] if pd.notna(row["fold_vs_null_mean_rate"]) else 0.0
        plt.annotate(str(row["feature"]), (row["rank"], y), fontsize=8, alpha=0.8)
    plt.axhline(1.0, linestyle="--", alpha=0.7)
    plt.title("Matched-background regulatory enrichment, fold versus null mean")
    plt.xlabel("Feature rank")
    plt.ylabel("Observed / null mean rate")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def write_logbook_entry(summary_json: dict, outpath: Path) -> None:
    strongest = summary_json["strongest_feature_by_p"]
    strongest_text = strongest if strongest is not None else "none"

    tex = rf"""
\subsection*{{Logbook Entry {ENTRY_NUMBER}: Matched-background regulatory enrichment test for the developmental bridge probes}}

\textbf{{What was being measured.}}
A matched-background enrichment test was performed to determine whether regulatory annotations observed among the focused developmental bridge probes were unusually concentrated relative to matched background probes from the same D3 versus not-D3 methylation universe.

\textbf{{Methods.}}
The focused probe set consisted of probes linked to PAX6, GLI2, BMP4, PTCH1, SHH, and OLIG2, extracted from the full-universe candidate evidence table. Background probes were drawn from the same D3 versus not-D3 universe after excluding probes belonging to the focused bridge genes and the local SIM2 anchor. Each focused probe was matched to background probes by absolute methylation-effect neighborhood and probe-context features, including refgene-group and CpG-island relation where possible. A total of {summary_json["n_matched_sets"]} matched background sets were generated. Enrichment was then tested for DNase labels, TFBS labels, open-chromatin labels, PHANTOM4 enhancer labels, enhancer-like regulatory names, and the \texttt{{Unclassified\_Cell\_type\_specific}} regulatory feature group.

\textbf{{Results.}}
The focused probe set contained {summary_json["n_focused_probes"]} probes. The strongest matched-background feature by empirical enrichment significance was: {latex_escape(strongest_text)}. Interpretation of this step depends on whether the observed regulatory annotations exceeded matched null expectations, rather than merely appearing somewhere among the bridge genes.

\textbf{{Tables written.}}
\texttt{{results/gse240704/matched\_background\_regulatory\_enrichment/matched\_background\_enrichment\_summary.csv}}\\
\texttt{{results/gse240704/matched\_background\_regulatory\_enrichment/matched\_background\_null\_summaries.csv}}\\
\texttt{{results/gse240704/matched\_background\_regulatory\_enrichment/tables\_tex/table\_matched\_background\_regulatory\_enrichment.tex}}\\
\texttt{{results/gse240704/matched\_background\_regulatory\_enrichment/matched\_background\_enrichment\_summary.json}}

\textbf{{Figures.}}
\texttt{{results/gse240704/matched\_background\_regulatory\_enrichment/plots/matched\_background\_feature\_enrichment.png}}

\textbf{{Interpretation.}}
This step determines whether the broader developmental bridge around the SIM2-centered signal has unusually concentrated regulatory annotation support, or whether the observed regulatory labels are compatible with matched background expectation. A positive enrichment result would strengthen the case for coordinated regulatory structure. A null result would favor interpretation as a distributed developmental backdrop without strong shared regulatory coordination in the current annotation layer.

\textbf{{Next steps.}}
Use the matched-background result to decide whether manuscript language should lean toward coordinated regulatory enrichment or remain conservative and describe only a weak distributed developmental backdrop around the dominant SIM2-local methylation event.
""".strip() + "\n"

    outpath.write_text(tex, encoding="utf-8")


def main() -> None:
    ensure_dir(OUTDIR)
    ensure_dir(PLOTDIR)
    ensure_dir(TABLETEXDIR)
    ensure_dir(LOGBOOKDIR)

    rng = random.Random(SEED)

    universe = load_universe_with_manifest()
    universe = add_matching_features(universe)

    focused_probe_ids = set(load_focused_probe_ids())
    cols = choose_cols(universe)
    probe_col = cols["probe_id"]
    if probe_col is None:
        raise ValueError("Could not identify probe_id column in universe table.")

    focused_df = universe[universe[probe_col].astype(str).isin(focused_probe_ids)].copy().reset_index(drop=True)
    if focused_df.empty:
        raise ValueError("Focused probe set is empty after matching to universe table.")

    background_df = build_background_pool(universe, focused_probe_ids)
    if background_df.empty:
        raise ValueError("Background pool is empty.")

    observed = summarize_feature_counts(focused_df)

    null_summaries = []
    for _ in range(N_MATCHED_SETS):
        sampled = sample_one_matched_set(focused_df, background_df, rng)
        null_summaries.append(summarize_feature_counts(sampled))

    enrichment_df = build_enrichment_summary(observed, null_summaries, n_probes=len(focused_df))
    enrichment_csv = OUTDIR / "matched_background_enrichment_summary.csv"
    enrichment_df.to_csv(enrichment_csv, index=False)

    null_df = pd.DataFrame(null_summaries)
    null_csv = OUTDIR / "matched_background_null_summaries.csv"
    null_df.to_csv(null_csv, index=False)

    tex_path = TABLETEXDIR / "table_matched_background_regulatory_enrichment.tex"
    write_summary_tex(enrichment_df, tex_path)

    plot_path = PLOTDIR / "matched_background_feature_enrichment.png"
    write_feature_plot(enrichment_df, plot_path)

    strongest_feature = None
    if not enrichment_df.empty:
        strongest_feature = enrichment_df.iloc[0]["feature"]

    summary_json = {
        "dataset": DATASET,
        "full_universe_table": str(FULL_UNIVERSE_TABLE),
        "focused_hits_table": str(FOCUSED_HITS_TABLE),
        "n_focused_probes": int(len(focused_df)),
        "n_background_pool_probes": int(len(background_df)),
        "n_matched_sets": N_MATCHED_SETS,
        "seed": SEED,
        "strongest_feature_by_p": strongest_feature,
        "outputs": {
            "enrichment_summary_csv": str(enrichment_csv),
            "null_summaries_csv": str(null_csv),
            "summary_tex": str(tex_path),
            "plot": str(plot_path),
        },
    }

    summary_json_path = OUTDIR / "matched_background_enrichment_summary.json"
    with open(summary_json_path, "w", encoding="utf-8") as f:
        json.dump(summary_json, f, indent=2)

    logbook_path = LOGBOOKDIR / "entry_33_matched_background_regulatory_enrichment.tex"
    write_logbook_entry(summary_json, logbook_path)

    print(json.dumps(summary_json, indent=2))
    print(f"[ok] wrote logbook entry to {logbook_path}")


if __name__ == "__main__":
    main()
