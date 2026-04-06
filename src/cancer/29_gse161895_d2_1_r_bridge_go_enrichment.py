#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def run_gseapy_enrichr(genes: list[str], out_dir: Path, label: str):
    try:
        import gseapy as gp
    except Exception as e:
        return {"status": "gseapy_not_available", "error": str(e), "output": None}

    if len(genes) == 0:
        return {"status": "no_genes", "output": None}

    enr_dir = out_dir / f"gseapy_{label}"
    enr_dir.mkdir(parents=True, exist_ok=True)

    try:
        enr = gp.enrichr(
            gene_list=genes,
            gene_sets=[
                "GO_Biological_Process_2023",
                "GO_Molecular_Function_2023",
                "GO_Cellular_Component_2023",
            ],
            organism="human",
            outdir=str(enr_dir),
            cutoff=1.0,
        )
        if enr is None or enr.results is None:
            return {"status": "ran_but_no_results", "output": None}

        res = enr.results.copy()
        out_csv = out_dir / f"{label}_go_enrichment_results.csv"
        res.to_csv(out_csv, index=False)
        return {"status": "ok", "output": str(out_csv)}
    except Exception as e:
        return {"status": "enrichr_failed", "error": str(e), "output": None}


def main() -> None:
    cfg = read_config()
    base_dir = REPO_ROOT / cfg["results_dir"] / "d2_1_r_treatment_bridge"
    out_dir = REPO_ROOT / cfg["results_dir"] / "d2_1_r_go_enrichment"
    out_dir.mkdir(parents=True, exist_ok=True)

    up = pd.read_csv(base_dir / "d2_1_r_treatment_bridge_up_hits.csv")
    down = pd.read_csv(base_dir / "d2_1_r_treatment_bridge_down_hits.csv")

    up_genes = sorted(set(up["gene_symbol"].dropna().astype(str)) - {"", "nan"})
    down_genes = sorted(set(down["gene_symbol"].dropna().astype(str)) - {"", "nan"})

    pd.DataFrame({"gene_symbol": up_genes}).to_csv(out_dir / "up_gene_list.csv", index=False)
    pd.DataFrame({"gene_symbol": down_genes}).to_csv(out_dir / "down_gene_list.csv", index=False)

    up_status = run_gseapy_enrichr(up_genes, out_dir, "up")
    down_status = run_gseapy_enrichr(down_genes, out_dir, "down")

    summary = {
        "up_gene_count": len(up_genes),
        "down_gene_count": len(down_genes),
        "up_status": up_status,
        "down_status": down_status,
    }

    with open(out_dir / "d2_1_r_go_enrichment_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
