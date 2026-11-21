import pandas as pd
import subprocess
import os

# 1. Helper: Translation function
def translate(codon):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    return table.get(codon.upper(), 'X')

# 2. Parse the annotated hotspots robustly
data = []
with open('MPXV_hotspots_annotated.txt', 'r') as f:
    for line in f:
        parts = line.split()
        if len(parts) < 6: continue
        # Chrom(0), Start(1), End(2), GeneID(5)
        # Note: GeneID usually looks like "ID=PHDGHHDA_00001"
        clean_id = parts[5].replace("ID=", "")
        data.append([parts[0], int(parts[1]), int(parts[2]), clean_id])

hotspots = pd.DataFrame(data, columns=['chrom', 'start', 'end', 'gene_id'])
unique_genes = hotspots['gene_id'].unique()

results = []
print(f"{'Gene_ID':<20} {'dN':<5} {'dS':<5} {'Ratio':<8} {'Selection'}")
print("-" * 65)

# Path to your prokka nucleotide file (adjust if name is different)
ffn_file = "ref/prokka_out/PROKKA_02232026.ffn" 

for gene in unique_genes:
    # A. Extract the DNA sequence for THIS gene using samtools
    # This ensures we are always translating the correct reading frame
    cmd_fa = f"samtools faidx {ffn_file} {gene}"
    res_fa = subprocess.run(cmd_fa, shell=True, capture_output=True, text=True)
    
    if res_fa.returncode != 0: continue
    # Strip the header and join the sequence
    ref_dna = "".join(res_fa.stdout.splitlines()[1:])
    
    # B. Get coordinates from our hotspot table
    gene_info = hotspots[hotspots['gene_id'] == gene].iloc[0]
    g_start = gene_info['start']
    
    # C. Extract and Load SNPs
    with open('current_gene.bed', 'w') as f:
        f.write(f"{gene_info['chrom']}\t{gene_info['start']}\t{gene_info['end']}\n")
    
    subprocess.run(f"bedtools intersect -a mpox_fixed.vcf -b current_gene.bed -wa > current_snps.vcf", shell=True)
    
    if not os.path.exists('current_snps.vcf') or os.path.getsize('current_snps.vcf') == 0:
        continue

    snps = pd.read_csv('current_snps.vcf', sep='\t', comment='#', header=None, dtype={1: int, 3: str, 4: str})
    snps = snps.drop_duplicates(subset=[1])

    # D. The Core Calculation
    dn, ds = 0, 0
    for _, row in snps.iterrows():
        abs_pos = row[1]
        alt_base = row[4]
        rel_pos = abs_pos - g_start - 1 # Map VCF to DNA string
        
        if rel_pos < 0 or rel_pos >= len(ref_dna): continue
        
        codon_start = (rel_pos // 3) * 3
        pos_in_codon = rel_pos % 3
        
        ref_codon = list(ref_dna[codon_start : codon_start+3])
        if len(ref_codon) < 3: continue
        
        alt_codon = ref_codon.copy()
        alt_codon[pos_in_codon] = alt_base
        
        ref_aa = translate("".join(ref_codon))
        alt_aa = translate("".join(alt_codon))
        
        if ref_aa == alt_aa: ds += 1
        else: dn += 1

    ratio = dn / ds if ds > 0 else (dn if dn > 0 else 0.0)
    sel = "Positive" if ratio > 1 else "Purifying"
    
    print(f"{gene:<20} {dn:<5} {ds:<5} {ratio:<8.2f} {sel}")
    results.append([gene, dn, ds, ratio, sel])

# 3. Final Output
df_final = pd.DataFrame(results, columns=['GeneID', 'dN', 'dS', 'Ratio', 'Selection'])
df_final.to_csv('all_hotspots_dnds.csv', index=False)