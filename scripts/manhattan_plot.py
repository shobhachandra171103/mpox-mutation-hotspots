import pandas as pd
import matplotlib.pyplot as plt
import re

# 1. Load the density data
df = pd.read_csv('MPXV_window_SNPdensity_500.tab', sep='\t', header=None,
                names=['chrom', 'start', 'end', 'count', 'c2', 'size', 'density'])
df['midpoint'] = (df['start'] + df['end']) / 2

# 2. Parse the Annotation file for labels
# We'll take the top 3 windows by count
annotations = []
try:
    with open('MPXV_hotspots_annotated.txt', 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 6:
                # Store: Start, Count, Gene_ID
                annotations.append({'start': int(parts[1]), 'count': int(parts[3]), 'label': parts[5]})
except FileNotFoundError:
    print("Annotation file not found, skipping labels.")

# Sort annotations by count to find the highest peaks
top_annotations = sorted(annotations, key=lambda x: x['count'], reverse=True)[:3]

# 3. Create the Plot
plt.figure(figsize=(15, 7))
scatter = plt.scatter(df['midpoint'], df['count'], 
                      c=df['count'], cmap='magma', 
                      s=30, alpha=0.8, edgecolor='none')

# 4. Add the labels for the top 3 peaks
for ann in top_annotations:
    plt.annotate(ann['label'], 
                 xy=(ann['start'] + 250, ann['count']), 
                 xytext=(0, 10), 
                 textcoords='offset points', 
                 ha='center', va='bottom',
                 fontsize=10, fontweight='bold',
                 arrowprops=dict(arrowstyle='->', color='black', lw=1))

# 5. Styling
plt.title('Annotated Manhattan Plot: MPXV Mutation Hotspots', fontsize=16)
plt.xlabel('Genome Position (bp)', fontsize=12)
plt.ylabel('SNP Count (500bp Window)', fontsize=12)
plt.axhline(y=5, color='gray', linestyle='--', alpha=0.3) # Baseline threshold
plt.colorbar(scatter, label='Mutation Intensity')
plt.tight_layout()

# 6. Save
plt.savefig('mpxv_annotated_manhattan.png', dpi=300)
print("Annotated Manhattan plot saved as mpxv_annotated_manhattan.png")
