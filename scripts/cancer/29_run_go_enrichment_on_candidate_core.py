import argparse
import json
from pathlib import Path

import requests
import yaml

from src.common.io import read_table, write_table


ENRICHR_ADD_URL = "https://maayanlab.cloud/Enrichr/addList"
ENRICHR_ENRICH_URL = "https://maayanlab.cloud/Enrichr/enrich"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--min-lists", type=int, default=3)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--library", default="GO_Biological_Process_2023")
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    tables_dir = Path(cfg["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    core = read_table(tables_dir / "irreversible_core_candidate_top_ranked.csv")
    core_sub = core[core["n_lists_present"] >= args.min_lists].head(args.top_n).copy()
    genes = core_sub["gene"].astype(str).tolist()

    if len(genes) == 0:
        raise ValueError("No candidate core genes selected for GO enrichment.")

    payload = {
        "list": (None, "\n".join(genes)),
        "description": (None, f"{cfg['dataset']['dataset_id']}_candidate_core"),
    }

    add_resp = requests.post(ENRICHR_ADD_URL, files=payload, timeout=120)
    add_resp.raise_for_status()
    add_data = add_resp.json()

    user_list_id = add_data["userListId"]

    enrich_resp = requests.get(
        ENRICHR_ENRICH_URL,
        params={
            "userListId": user_list_id,
            "backgroundType": args.library,
        },
        timeout=120,
    )
    enrich_resp.raise_for_status()
    enrich_data = enrich_resp.json()

    rows = []
    for entry in enrich_data.get(args.library, []):
        # Enrichr format:
        # [rank, term_name, p_value, z_score, combined_score, overlapping_genes, adjusted_p_value, old_p_value, old_adjusted_p_value]
        rows.append(
            {
                "rank": entry[0],
                "term_name": entry[1],
                "p_value": entry[2],
                "z_score": entry[3],
                "combined_score": entry[4],
                "overlapping_genes": ";".join(entry[5]) if isinstance(entry[5], list) else str(entry[5]),
                "adjusted_p_value": entry[6],
                "old_p_value": entry[7] if len(entry) > 7 else None,
                "old_adjusted_p_value": entry[8] if len(entry) > 8 else None,
            }
        )

    import pandas as pd
    res = pd.DataFrame(rows).sort_values(["adjusted_p_value", "p_value", "combined_score"], ascending=[True, True, False])

    out_path = tables_dir / "candidate_core_go_enrichment.csv"
    write_table(res, out_path)

    summary = {
        "dataset_id": cfg["dataset"]["dataset_id"],
        "library": args.library,
        "min_lists": int(args.min_lists),
        "top_n": int(args.top_n),
        "n_input_genes": int(len(genes)),
        "input_genes": genes,
        "user_list_id": int(user_list_id),
        "output_file": str(out_path),
        "n_terms": int(res.shape[0]),
    }

    with open(tables_dir / "candidate_core_go_enrichment_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
