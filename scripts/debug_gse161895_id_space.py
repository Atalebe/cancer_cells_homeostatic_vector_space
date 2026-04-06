#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import json

REPO = Path.home() / "cancer_cells_hrsm_sims"
DATASET = "gse161895"

state = pd.read_parquet(REPO / "data" / "processed" / DATASET / "state_table.parquet")
meta = pd.read_parquet(REPO / "data" / "metadata" / DATASET / "cell_metadata_registry.parquet")

print("\n[state columns]")
print(state.columns.tolist())

print("\n[meta columns]")
print(meta.columns.tolist())

print("\n[state head]")
print(state.head(10).to_string())

print("\n[meta head selected]")
cols = [c for c in ["cell_id", "gsm", "geo_accession", "title", "patient_or_donor_id", "source_class", "treatment_state"] if c in meta.columns]
print(meta[cols].head(10).to_string())

if "cell_id" in state.columns:
    s = state["cell_id"].astype(str)
    print("\n[first 30 state cell_id values]")
    print(json.dumps(s.head(30).tolist(), indent=2))

for c in ["cell_id", "gsm", "geo_accession", "title"]:
    if c in meta.columns:
        print(f"\n[first 30 meta {c} values]")
        print(json.dumps(meta[c].astype(str).head(30).tolist(), indent=2))
