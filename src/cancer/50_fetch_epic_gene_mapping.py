import pandas as pd
from pybiomart import Dataset

# 1. Load your top 200 probes
probes_df = pd.read_csv('results/gse240704/s_axis_drivers/s_axis_top_probe_drivers_annotated.csv')
probe_list = probes_df['ID_REF'].unique().tolist()

print(f"[*] Querying BioMart for {len(probe_list)} probes...")

# 2. Connect to Ensembl Human Dataset (GRCh38)
dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')

# 3. Run Query
# Note: 'illumina_methylation_epic' is the filter for EPIC/850k probes
results = dataset.query(
    attributes=['illumina_methylation_epic', 'external_gene_name', 'chromosome_name', 'start_position', 'transcript_biotype'],
    filters={'illumina_methylation_epic': probe_list}
)

# 4. Cleanup and Save
results.columns = ['probe_id', 'gene_symbol', 'chr', 'pos', 'type']
results = results.dropna(subset=['gene_symbol'])

out_path = 'results/gse240704/s_axis_drivers/s_axis_drivers_gene_mapped.csv'
results.to_csv(out_path, index=False)

print(f"[ok] Mapped {results['probe_id'].nunique()} probes to {results['gene_symbol'].nunique()} genes.")
print(f"[info] Results saved to: {out_path}")
