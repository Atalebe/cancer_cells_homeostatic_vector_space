import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

from src.common.io import read_table, write_table


def load_ranked_gene_list(path: Path, gene_col: str = "gene") -> pd.DataFrame:
    df = read_table(path).copy()
    if gene_col not in df.columns:
        raise ValueError(f"Missing gene column in {path}")
    df[gene_col] = df[gene_col].astype(str)
    df = df.reset_index(drop=True)
    df["rank"] = df.index + 1
    return df[[gene_col, "rank"]]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--top-n", type=int, default=50)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    tables_dir = Path(cfg["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    input_files = {
        "positive_H": tables_dir / "top_positive_genes_corr_H.csv",
        "negative_R": tables_dir / "top_negative_genes_corr_R.csv",
        "malignant_reservoir": tables_dir / "top_up_genes_malignant_reservoir.csv",
        "unstable_committed_malignant": tables_dir / "top_up_genes_unstable_committed_malignant.csv",
        "PKH26_High_vs_G1": tables_dir / "top_up_genes_PKH26_High_vs_G1.csv",
    }

    missing = [name for name, path in input_files.items() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required input ranking tables: {missing}")

    rank_maps = {}
    gene_sets = {}

    for name, path in input_files.items():
        df = load_ranked_gene_list(path)
        df = df.head(args.top_n).copy()
        rank_maps[name] = dict(zip(df["gene"], df["rank"]))
        gene_sets[name] = set(df["gene"])

    all_genes = sorted(set().union(*gene_sets.values()))

    rows = []
    for gene in all_genes:
        present_in = [name for name, genes in gene_sets.items() if gene in genes]
        n_lists = len(present_in)

        row = {
            "gene": gene,
            "n_lists_present": n_lists,
            "present_in_lists": ";".join(present_in),
        }

        ranks = []
        for name in input_files:
            rank_val = rank_maps[name].get(gene)
            row[f"rank_{name}"] = rank_val if rank_val is not None else ""
            if rank_val is not None:
                ranks.append(rank_val)

        row["mean_rank_across_present_lists"] = float(sum(ranks) / len(ranks)) if ranks else None
        row["priority_score"] = float(n_lists * 1000 - (sum(ranks) / len(ranks))) if ranks else 0.0
        rows.append(row)

    out_df = pd.DataFrame(rows).sort_values(
        ["n_lists_present", "mean_rank_across_present_lists"],
        ascending=[False, True]
    ).reset_index(drop=True)

    out_path = tables_dir / "irreversible_core_candidate_overlap.csv"
    write_table(out_df, out_path)

    # subsets by recurrence threshold
    threshold_outputs = {}
    for k in [2, 3, 4, 5]:
        sub = out_df[out_df["n_lists_present"] >= k].copy()
        sub_path = tables_dir / f"irreversible_core_candidates_at_least_{k}_lists.csv"
        write_table(sub, sub_path)
        threshold_outputs[f"at_least_{k}_lists"] = str(sub_path)

    top_ranked = out_df.copy().sort_values(
        ["n_lists_present", "priority_score"],
        ascending=[False, False]
    ).reset_index(drop=True)

    top_ranked_path = tables_dir / "irreversible_core_candidate_top_ranked.csv"
    write_table(top_ranked, top_ranked_path)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "top_n_from_each_input_list": int(args.top_n),
        "input_files": {k: str(v) for k, v in input_files.items()},
        "outputs": {
            "overlap_table": str(out_path),
            "top_ranked": str(top_ranked_path),
            "threshold_tables": threshold_outputs,
        },
        "n_unique_candidate_genes": int(out_df.shape[0]),
        "n_candidates_at_least_2_lists": int((out_df["n_lists_present"] >= 2).sum()),
        "n_candidates_at_least_3_lists": int((out_df["n_lists_present"] >= 3).sum()),
        "n_candidates_at_least_4_lists": int((out_df["n_lists_present"] >= 4).sum()),
        "n_candidates_all_5_lists": int((out_df["n_lists_present"] >= 5).sum()),
    }

    with open(tables_dir / "irreversible_core_candidate_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
