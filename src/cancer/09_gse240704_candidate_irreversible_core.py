from __future__ import annotations

from pathlib import Path

import pandas as pd

from _gse240704_utils import load_config


def main(config_path: str) -> None:
    cfg = load_config(config_path)

    state_path = cfg["outputs"]["state_table"]
    domain_path = "results/gse240704/state_domains/state_domains.parquet"
    out_dir = Path("results/gse240704/irreversible_core")
    out_dir.mkdir(parents=True, exist_ok=True)

    state = pd.read_parquet(state_path)

    if Path(domain_path).exists():
        dom = pd.read_parquet(domain_path)[["sample_id", "state_domain"]]
        state = state.merge(dom, on="sample_id", how="left")

    state["phi_q90_flag"] = state["phi"] >= state["phi"].quantile(0.90)
    state["phi_q95_flag"] = state["phi"] >= state["phi"].quantile(0.95)
    state["S_high_flag"] = state["S"] >= state["S"].quantile(0.90)
    state["M_high_flag"] = state["M"] >= state["M"].quantile(0.90)
    state["R_low_flag"] = state["R"] <= state["R"].quantile(0.10)

    state["candidate_stable_shell"] = (
        state["phi_q90_flag"] &
        state["S_high_flag"] &
        state["M_high_flag"]
    )

    state["candidate_irreversible_shell"] = (
        state["phi_q95_flag"] &
        state["S_high_flag"] &
        state["M_high_flag"] &
        state["R_low_flag"]
    )

    ranked = state.sort_values(
        ["candidate_irreversible_shell", "candidate_stable_shell", "phi", "S", "M"],
        ascending=[False, False, False, False, False]
    ).copy()

    ranked.to_parquet(out_dir / "candidate_irreversible_shell_table.parquet", index=False)
    ranked.to_csv(out_dir / "candidate_irreversible_shell_table.csv", index=False)

    summary = pd.DataFrame(
        [
            {"metric": "n_samples", "value": int(len(ranked))},
            {"metric": "candidate_stable_shell_n", "value": int(ranked["candidate_stable_shell"].sum())},
            {"metric": "candidate_irreversible_shell_n", "value": int(ranked["candidate_irreversible_shell"].sum())},
            {"metric": "phi_q90", "value": float(state["phi"].quantile(0.90))},
            {"metric": "phi_q95", "value": float(state["phi"].quantile(0.95))},
            {"metric": "S_q90", "value": float(state["S"].quantile(0.90))},
            {"metric": "M_q90", "value": float(state["M"].quantile(0.90))},
            {"metric": "R_q10", "value": float(state["R"].quantile(0.10))},
        ]
    )
    summary.to_csv(out_dir / "candidate_irreversible_shell_summary.csv", index=False)

    print(f"[ok] wrote candidate shell tables to {out_dir}")
    print(f"[info] candidate stable shell n: {int(ranked['candidate_stable_shell'].sum())}")
    print(f"[info] candidate irreversible shell n: {int(ranked['candidate_irreversible_shell'].sum())}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    args = parser.parse_args()
    main(args.config)
