import pandas as pd
import matplotlib.pyplot as plt
import re

data = []
# Manually parse to handle the space-heavy 'product' column
with open('MPXV_hotspots_annotated.txt', 'r') as f:
    for line in f:
        parts = line.split()
        if len(parts) < 6: continue
        # Extract fixed columns: Chrom, Start, End, Count, Density, ID
        # Everything after index 5 is the 'Product' description
        chrom, start, end, count, density, gene_id = parts[0:6]
        product = " ".join(parts[6:])
        data.append([chrom, int(start), int(end), int(count), float(density.replace('Density:', '')), gene_id, product])

df = pd.DataFrame(data, columns=['chrom', 'start', 'end', 'count', 'density', 'gene_id', 'product'])

# Group by Gene ID and Product to get the total SNP count per gene
gene_hits = df.groupby(['gene_id', 'product'])['count'].sum().sort_values().tail(10)

# Create the plot
plt.figure(figsize=(12, 8))
y_labels = [f"{idx[0]} ({idx[1]})" for idx in gene_hits.index]
plt.hlines(y=y_labels, xmin=0, xmax=gene_hits.values, color='firebrick', linewidth=2)
plt.plot(gene_hits.values, y_labels, "o", color='firebrick', markersize=10)

plt.title('Top 10 Most Mutated Genes (Total SNPs in Hotspot Windows)', fontsize=14)
plt.xlabel('Cumulative SNP Count', fontsize=12)
plt.ylabel('Gene ID and Product', fontsize=12)
plt.grid(axis='x', linestyle='--', alpha=0.6)
plt.tight_layout()

plt.savefig('mpxv_top_genes_fixed.png')
print("Success! Graph saved as mpxv_top_genes_fixed.png")
