from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from common.score_utils import robust_zscore, weighted_sum


@dataclass
class ProxyRecipe:
    h_weights: dict[str, float]
    s_weights: dict[str, float]
    m_weights: dict[str, float]
    r_weights: dict[str, float]


def build_state_coordinates(feature_table, recipe: ProxyRecipe):
    out = feature_table.copy()
    out["H"] = robust_zscore(weighted_sum(out, recipe.h_weights))
    out["S"] = robust_zscore(weighted_sum(out, recipe.s_weights))
    out["M"] = robust_zscore(weighted_sum(out, recipe.m_weights))
    out["R"] = robust_zscore(weighted_sum(out, recipe.r_weights))
    return out
