import os
import GEOparse
import pandas as pd

# Setup paths
GSE_ID = "GSE240704"
DATA_DIR = "data/raw/annotations"
os.makedirs(DATA_DIR, exist_ok=True)

print(f"[*] Fetching platform metadata for {GSE_ID}...")

try:
    # 1. Download GSE to identify the platform (GPL)
    gse = GEOparse.get_GEO(geo=GSE_ID, destdir=DATA_DIR, silent=True)
    
    # Identify the primary platform ID
    gpl_id = list(gse.gpls.keys())[0]
    print(f"[+] Identified Platform: {gpl_id}")

    # 2. Access the GPL SOFT file table
    gpl = gse.gpls[gpl_id]
    annotation_df = gpl.table # Contains ID, Gene Symbol, etc.
    
    # 3. Clean and Save
    # Keep only essential columns to save space if needed
    cols_to_keep = ['ID', 'UCSC_RefGene_Name', 'UCSC_RefGene_Group', 'UCSC_CpG_Islands_Name', 'Chromosome', 'Start']
    # Filtering for columns that actually exist in this specific GPL
    existing_cols = [c for c in cols_to_keep if c in annotation_df.columns]
    
    out_path = os.path.join(DATA_DIR, f"{gpl_id}_annotation.parquet")
    annotation_df[existing_cols].to_parquet(out_path)
    
    print(f"[ok] Annotation saved to: {out_path}")
    print(f"[info] Columns mapped: {existing_cols}")
    print(f"[info] Total probes: {len(annotation_df)}")

except Exception as e:
    print(f"[error] Failed to fetch: {e}")
