#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import pandas as pd


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
DATASET = "gse240704"

INPUT_CSV = (
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
    / "sim2_followup"
    / "tables_tex"
)

OUT_TEX = OUTDIR / "table_sim2_downstream_target_seed.tex"


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


def main() -> None:
    ensure_dir(OUTDIR)

    if not INPUT_CSV.exists():
        raise FileNotFoundError(f"Missing input csv: {INPUT_CSV}")

    df = pd.read_csv(INPUT_CSV)

    required = ["gene", "category", "note"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Curated seed list of candidate downstream or pathway-linked genes used for SIM2-centered interpretive follow up. This table is a hypothesis scaffold for later mechanistic overlay and does not itself establish direct regulation in GSE240704.}")
    lines.append(r"\label{tab:sim2_downstream_target_seed}")
    lines.append(r"\begin{tabular}{lll}")
    lines.append(r"\hline")
    lines.append(r"Gene & Category & Note \\")
    lines.append(r"\hline")

    for _, row in df.iterrows():
        gene = latex_escape(row["gene"])
        category = latex_escape(row["category"])
        note = latex_escape(row["note"])
        lines.append(f"{gene} & {category} & {note} \\\\")

    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    OUT_TEX.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[ok] wrote {OUT_TEX}")


if __name__ == "__main__":
    main()
