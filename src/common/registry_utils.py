from __future__ import annotations

from pathlib import Path
import pandas as pd

DATASET_REGISTRY_COLUMNS = [
    "dataset_id",
    "source_type",
    "accession",
    "title",
    "tumor_type",
    "species",
    "assay",
    "branch_target",
    "has_processed_matrix",
    "has_metadata",
    "has_cell_annotations",
    "has_stemness_annotation",
    "has_treatment_annotation",
    "has_normal_reference",
    "has_timecourse",
    "has_resistance_label",
    "n_cells_reported",
    "n_genes_reported",
    "priority_tier",
    "expected_H_proxy",
    "expected_S_proxy",
    "expected_M_proxy",
    "expected_R_proxy",
    "cell_cycle_confound_risk",
    "metabolic_confound_risk",
    "batch_risk",
    "duplication_risk",
    "status",
    "local_raw_dir",
    "local_processed_dir",
    "notes",
]


def empty_dataset_registry() -> pd.DataFrame:
    return pd.DataFrame(columns=DATASET_REGISTRY_COLUMNS)


def write_registry(path: str | Path) -> None:
    df = empty_dataset_registry()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
