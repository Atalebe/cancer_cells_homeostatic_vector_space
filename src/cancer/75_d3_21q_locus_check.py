#!/usr/bin/env python3

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd


@dataclass
class Locus:
    gene: str
    chromosome: str
    start: int
    end: int


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

CENTROMERE_TABLE = (
    REPO_ROOT
    / "data"
    / "reference"
    / "cytoband"
    / "hg19_centromere_breakpoints.csv"
)

OUTDIR = REPO_ROOT / "results" / DATASET / "chr21_locus_followup"

WINDOW_BP = 2_000_000

# hg19 approximate locus ranges, good enough for local inspection windows
TARGET_LOCI = [
    Locus(gene="ERG", chromosome="21", start=39_800_000, end=40_300_000),
    Locus(gene="RUNX1", chromosome="21", start=36_150_000, end=36_450_000),
]


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
    """
    Keep only the first instance of duplicated column names.
    This avoids pandas returning a DataFrame when a repeated name is selected.
    """
    return df.loc[:, ~df.columns.duplicated()].copy()


def find_first_existing(df: pd.DataFrame, candidates: list[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def get_series(df: pd.DataFrame, colname: str) -> pd.Series:
    x = df.loc[:, colname]
    if isinstance(x, pd.DataFrame):
        return x.iloc[:, 0]
    return x


def find_probe_id_column(df: pd.DataFrame) -> str:
    candidates = [
        "probe_id",
        "ilmnid",
        "ilmn_id",
        "name",
        "cg_id",
        "probe",
        "cpg_id",
        "id_ref",
        "id",
    ]
    col = find_first_existing(df, candidates)
    if col is None:
        raise ValueError(f"Could not find probe id column. Columns: {list(df.columns)}")
    return col


def find_direction_column(df: pd.DataFrame) -> str:
    candidates = ["direction", "probe_direction", "contrast_direction"]
    col = find_first_existing(df, candidates)
    if col is None:
        raise ValueError(f"Could not find direction column. Columns: {list(df.columns)}")
    return col


def choose_chr_column(df: pd.DataFrame) -> str:
    candidates = [
        "chromosome_ann",
        "chromosome_dir",
        "chromosome",
        "chr_ann",
        "chr_dir",
        "chr",
        "chrom_ann",
        "chrom_dir",
        "chrom",
    ]
    col = find_first_existing(df, candidates)
    if col is None:
        raise ValueError(f"Could not find chromosome column. Columns: {list(df.columns)}")
    return col


def choose_pos_column(df: pd.DataFrame) -> str:
    candidates = [
        "mapinfo_ann",
        "mapinfo_dir",
        "mapinfo",
        "position_ann",
        "position_dir",
        "position",
        "pos_ann",
        "pos_dir",
        "pos",
        "bp_ann",
        "bp_dir",
        "bp",
    ]
    col = find_first_existing(df, candidates)
    if col is None:
        raise ValueError(f"Could not find genomic position column. Columns: {list(df.columns)}")
    return col


def choose_gene_column(df: pd.DataFrame) -> Optional[str]:
    candidates = [
        "ucsc_refgene_name_ann",
        "ucsc_refgene_name_dir",
        "ucsc_refgene_name",
        "gene_symbol",
        "gene_symbols",
        "genes",
        "refgene_name",
    ]
    return find_first_existing(df, candidates)


def choose_refgene_group_column(df: pd.DataFrame) -> Optional[str]:
    candidates = [
        "refgene_group_ann",
        "refgene_group_dir",
        "ucsc_refgene_group_ann",
        "ucsc_refgene_group_dir",
        "refgene_group",
        "ucsc_refgene_group",
    ]
    return find_first_existing(df, candidates)


def choose_island_column(df: pd.DataFrame) -> Optional[str]:
    candidates = [
        "relation_to_cpg_island_ann",
        "relation_to_cpg_island_dir",
        "relation_to_ucsc_cpg_island_ann",
        "relation_to_ucsc_cpg_island_dir",
        "relation_to_cpg_island",
        "relation_to_ucsc_cpg_island",
    ]
    return find_first_existing(df, candidates)


def choose_effect_column(df: pd.DataFrame) -> Optional[str]:
    candidates = [
        "delta_median_a_minus_b",
        "delta_mean_a_minus_b",
        "abs_delta_median",
        "rank_biserial",
        "abs_rank_biserial",
        "effect_size",
        "loading",
        "score",
        "t_stat",
        "logfc",
    ]
    return find_first_existing(df, candidates)


def standardize_chr_value(x: object) -> Optional[str]:
    if pd.isna(x):
        return None

    s = str(x).strip()
    if not s:
        return None

    s = re.sub(r"^chr", "", s, flags=re.IGNORECASE).strip()

    if s.upper() in {"X", "Y"}:
        return s.upper()

    if s.upper() in {"M", "MT"}:
        return "MT"

    if s.isdigit():
        return str(int(s))

    return s


def split_gene_tokens(value: object) -> list[str]:
    if pd.isna(value):
        return []

    s = str(value).strip()
    if not s:
        return []

    tokens: list[str] = []
    for part in re.split(r"[;,/|]+", s):
        token = part.strip()
        if token:
            tokens.append(token)
    return tokens


def load_centromere_breakpoints(path: Path) -> dict[str, float]:
    cen = pd.read_csv(path)
    cen = normalize_columns(cen)
    cen = deduplicate_columns(cen)

    chr_col = find_first_existing(cen, ["chromosome", "chr", "chrom"])
    if chr_col is None:
        raise ValueError(
            f"Centromere file missing chromosome column. Columns: {list(cen.columns)}"
        )

    bp_col = find_first_existing(cen, ["centromere_bp", "centromere", "q_start", "p_end"])
    if bp_col is None:
        raise ValueError(
            f"Centromere file missing breakpoint column. Columns: {list(cen.columns)}"
        )

    out = cen[[chr_col, bp_col]].copy()
    out.columns = ["chromosome", "centromere_bp"]
    out["chromosome"] = out["chromosome"].map(standardize_chr_value)
    out["centromere_bp"] = pd.to_numeric(out["centromere_bp"], errors="coerce")
    out = out.dropna(subset=["chromosome", "centromere_bp"]).drop_duplicates()

    return dict(zip(out["chromosome"], out["centromere_bp"]))


def assign_arm(chrom: object, pos: object, centromere_map: dict[str, float]) -> Optional[str]:
    chrom_std = standardize_chr_value(chrom)
    if chrom_std is None or pd.isna(pos):
        return None
    if chrom_std not in centromere_map:
        return None
    return f"{chrom_std}p" if float(pos) < float(centromere_map[chrom_std]) else f"{chrom_std}q"


def count_gene_tokens(series: pd.Series) -> pd.DataFrame:
    counts: dict[str, int] = {}

    for value in series.fillna("").astype(str):
        for token in split_gene_tokens(value):
            counts[token] = counts.get(token, 0) + 1

    if not counts:
        return pd.DataFrame(columns=["gene_token", "n_probes"])

    out = pd.DataFrame(
        {
            "gene_token": list(counts.keys()),
            "n_probes": list(counts.values()),
        }
    ).sort_values(["n_probes", "gene_token"], ascending=[False, True])

    return out.reset_index(drop=True)


def build_window_table(
    df: pd.DataFrame,
    locus: Locus,
    window_bp: int,
    gene_col: Optional[str],
) -> pd.DataFrame:
    start = locus.start - window_bp
    end = locus.end + window_bp

    out = df[
        (df["chromosome_std"] == locus.chromosome)
        & (df["position_std"].notna())
        & (df["position_std"] >= start)
        & (df["position_std"] <= end)
    ].copy()

    out["target_gene"] = locus.gene
    out["locus_start"] = locus.start
    out["locus_end"] = locus.end
    out["window_start"] = start
    out["window_end"] = end
    out["distance_to_locus_mid_bp"] = out["position_std"] - ((locus.start + locus.end) / 2.0)

    if gene_col is not None and gene_col in out.columns:
        gene_series = get_series(out, gene_col).fillna("").astype(str)
        out["contains_target_gene_token"] = gene_series.map(
            lambda s: locus.gene in split_gene_tokens(s)
        )
    else:
        out["contains_target_gene_token"] = False

    effect_col = choose_effect_column(out)

    if effect_col is not None and effect_col in out.columns:
        out["_abs_effect"] = pd.to_numeric(get_series(out, effect_col), errors="coerce").abs()
        out = out.sort_values(
            ["contains_target_gene_token", "_abs_effect", "position_std"],
            ascending=[False, False, True],
        )
    else:
        out = out.sort_values(
            ["contains_target_gene_token", "position_std"],
            ascending=[False, True],
        )

    return out.reset_index(drop=True)


def main() -> None:
    ensure_dir(OUTDIR)

    if not DIRECTIONAL_TABLE.exists():
        raise FileNotFoundError(f"Missing directional table: {DIRECTIONAL_TABLE}")

    if not MANIFEST_TABLE.exists():
        raise FileNotFoundError(f"Missing manifest table: {MANIFEST_TABLE}")

    if not CENTROMERE_TABLE.exists():
        raise FileNotFoundError(f"Missing centromere table: {CENTROMERE_TABLE}")

    directional = pd.read_csv(DIRECTIONAL_TABLE)
    manifest = pd.read_parquet(MANIFEST_TABLE)

    directional = normalize_columns(directional)
    manifest = normalize_columns(manifest)

    directional = deduplicate_columns(directional)
    manifest = deduplicate_columns(manifest)

    dir_probe_col = find_probe_id_column(directional)
    man_probe_col = find_probe_id_column(manifest)
    direction_col = find_direction_column(directional)

    directional[dir_probe_col] = directional[dir_probe_col].astype(str).str.strip()
    manifest[man_probe_col] = manifest[man_probe_col].astype(str).str.strip()

    merged = directional.merge(
        manifest,
        left_on=dir_probe_col,
        right_on=man_probe_col,
        how="left",
        suffixes=("_dir", "_ann"),
    )

    merged = normalize_columns(merged)
    merged = deduplicate_columns(merged)

    chr_col = choose_chr_column(merged)
    pos_col = choose_pos_column(merged)
    gene_col = choose_gene_column(merged)
    refgene_group_col = choose_refgene_group_column(merged)
    island_col = choose_island_column(merged)
    effect_col = choose_effect_column(merged)

    merged["chromosome_std"] = get_series(merged, chr_col).map(standardize_chr_value)
    merged["position_std"] = pd.to_numeric(get_series(merged, pos_col), errors="coerce")

    centromere_map = load_centromere_breakpoints(CENTROMERE_TABLE)

    merged["chromosome_arm"] = [
        assign_arm(ch, pos, centromere_map)
        for ch, pos in zip(merged["chromosome_std"], merged["position_std"])
    ]

    lower = merged[
        get_series(merged, direction_col).astype(str).str.lower().eq("lower_in_a")
    ].copy()

    chr21 = lower[lower["chromosome_std"] == "21"].copy()
    chr21q = chr21[chr21["chromosome_arm"] == "21q"].copy()

    if effect_col is not None and effect_col in chr21q.columns:
        chr21q["_abs_effect"] = pd.to_numeric(
            get_series(chr21q, effect_col), errors="coerce"
        ).abs()
        chr21q = chr21q.sort_values(
            ["_abs_effect", "position_std"],
            ascending=[False, True],
        ).reset_index(drop=True)
    else:
        chr21q = chr21q.sort_values("position_std").reset_index(drop=True)

    probe_table_path = OUTDIR / "D3_21q_probe_table.csv"
    chr21q.to_csv(probe_table_path, index=False)

    if gene_col is not None and gene_col in chr21q.columns:
        gene_counts_df = count_gene_tokens(get_series(chr21q, gene_col))
    else:
        gene_counts_df = pd.DataFrame(columns=["gene_token", "n_probes"])

    gene_counts_path = OUTDIR / "D3_21q_gene_token_counts.csv"
    gene_counts_df.to_csv(gene_counts_path, index=False)

    locus_rows = []
    for locus in TARGET_LOCI:
        win = build_window_table(
            chr21q,
            locus=locus,
            window_bp=WINDOW_BP,
            gene_col=gene_col,
        )

        out_csv = OUTDIR / f"D3_21q_{locus.gene}_probe_window.csv"
        win.to_csv(out_csv, index=False)

        locus_rows.append(
            {
                "gene": locus.gene,
                "chromosome": locus.chromosome,
                "locus_start": locus.start,
                "locus_end": locus.end,
                "window_bp_each_side": WINDOW_BP,
                "n_window_probes": int(len(win)),
                "n_direct_gene_token_hits": int(win["contains_target_gene_token"].sum())
                if "contains_target_gene_token" in win.columns
                else 0,
            }
        )

    locus_summary_df = pd.DataFrame(locus_rows)
    locus_summary_path = OUTDIR / "D3_21q_locus_compact_summary.csv"
    locus_summary_df.to_csv(locus_summary_path, index=False)

    context_summary: dict[str, dict] = {}

    if refgene_group_col is not None and refgene_group_col in chr21q.columns:
        context_summary["refgene_group_counts"] = (
            get_series(chr21q, refgene_group_col)
            .fillna("NA")
            .astype(str)
            .value_counts()
            .to_dict()
        )

    if island_col is not None and island_col in chr21q.columns:
        context_summary["cpg_island_relation_counts"] = (
            get_series(chr21q, island_col)
            .fillna("NA")
            .astype(str)
            .value_counts()
            .to_dict()
        )

    summary = {
        "dataset": DATASET,
        "direction_filter": "lower_in_a",
        "n_lower_in_a_total": int(len(lower)),
        "n_chr21_lower_in_a": int(len(chr21)),
        "n_chr21q_lower_in_a": int(len(chr21q)),
        "directional_table": str(DIRECTIONAL_TABLE),
        "manifest_table": str(MANIFEST_TABLE),
        "centromere_table": str(CENTROMERE_TABLE),
        "selected_columns": {
            "direction": direction_col,
            "chromosome": chr_col,
            "position": pos_col,
            "gene": gene_col,
            "refgene_group": refgene_group_col,
            "cpg_island_relation": island_col,
            "effect": effect_col,
        },
        "outputs": {
            "probe_table": str(probe_table_path),
            "gene_token_counts": str(gene_counts_path),
            "locus_compact_summary": str(locus_summary_path),
            "erg_window": str(OUTDIR / "D3_21q_ERG_probe_window.csv"),
            "runx1_window": str(OUTDIR / "D3_21q_RUNX1_probe_window.csv"),
        },
        "context_summary": context_summary,
        "locus_summary_rows": locus_rows,
    }

    summary_path = OUTDIR / "D3_21q_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
