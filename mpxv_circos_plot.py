import pandas as pd
from pycirclize import Circos
import matplotlib.pyplot as plt

# 1. Load and prepare data
dnds_df = pd.read_csv('all_hotspots_dnds.csv')

pos_data = []
with open('MPXV_hotspots_annotated.txt', 'r') as f:
    for line in f:
        parts = line.split()
        if len(parts) < 6: continue
        clean_id = parts[5].replace("ID=", "")
        pos_data.append([clean_id, int(parts[1]), int(parts[2])])

pos_df = pd.DataFrame(pos_data, columns=['GeneID', 'Start', 'End'])
df = pd.merge(dnds_df, pos_df, on='GeneID')

# 2. Initialize Circos
genome_length = int(df['End'].max() + 5000) 
circos = Circos(sectors={"MPXV": genome_length}, space=5)

for sector in circos.sectors:
    # --- OUTER TRACK: dN/dS Ratios (The Selection Pressure) ---
    res_track = sector.add_track((80, 100))
    res_track.axis()
    
    x_pos = df['Start'].tolist()
    ratios = [min(r, 10) for r in df['Ratio']] 
    colors = ["#e74c3c" if r > 1.0 else "#3498db" for r in df['Ratio']]
    
    res_track.bar(x_pos, ratios, width=1200, color=colors)
    res_track.text("dN/dS Ratio", r=110, size=12)

    # --- INNER TRACK: Mutation Count (The SNP Intensity) ---
    den_track = sector.add_track((55, 75))
    den_track.axis()
    
    max_muts = (df['dN'] + df['dS']).max()
    for _, row in df.iterrows():
        total = row['dN'] + row['dS']
        a = max(0.2, min(total / max_muts, 0.8)) if max_muts > 0 else 0.5
        den_track.rect(row['Start'], row['End'], color="orange", alpha=a)
    
    den_track.text("Gene Density", r=45, size=12)

    # --- ANNOTATIONS: Labeling the outliers ---
    # We use track.text for a more stable API approach
    top_hits = df.nlargest(5, 'Ratio')
    for _, row in top_hits.iterrows():
        label = f"{row['GeneID'].split('_')[-1]} ({row['Ratio']:.1f})"
        # This places text at a specific radius (r) and genomic position
        res_track.text(label, x=row['Start'], r=120, size=9, orientation="vertical")

# 3. Save the result
circos.savefig("mpxv_circular_analysis.png")
print(f"Success! Circular plot saved for MPXV genome ({genome_length} bp).")