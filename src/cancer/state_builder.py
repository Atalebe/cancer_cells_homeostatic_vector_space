from __future__ import annotations

from pathlib import Path
import yaml
import pandas as pd

from common.io import read_table, write_table
from cancer.proxy_definitions import ProxyRecipe, build_state_coordinates


def load_recipe_from_yaml(path: str | Path) -> ProxyRecipe:
    with open(path, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    return ProxyRecipe(
        h_weights=cfg["H"]["weights"],
        s_weights=cfg["S"]["weights"],
        m_weights=cfg["M"]["weights"],
        r_weights=cfg["R"]["weights"],
    )


def build_state_table(feature_table_path: str | Path, proxy_recipe_path: str | Path, output_path: str | Path) -> None:
    feature_table = read_table(feature_table_path)
    recipe = load_recipe_from_yaml(proxy_recipe_path)
    state_table = build_state_coordinates(feature_table, recipe)
    write_table(state_table, output_path)
