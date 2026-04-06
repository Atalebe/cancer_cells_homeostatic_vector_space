from __future__ import annotations

import pandas as pd


def simple_branch_summary(state_table: pd.DataFrame, group_col: str) -> pd.DataFrame:
    metrics = ["H", "S", "M", "R"]
    available = [m for m in metrics if m in state_table.columns]
    return state_table.groupby(group_col)[available].agg(["mean", "median", "count"])
