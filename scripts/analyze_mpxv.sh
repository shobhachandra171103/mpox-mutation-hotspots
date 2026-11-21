#!/bin/bash

# --- 1. SET FILENAMES ---
VCF="mpox.vcf"
FIXED_VCF="mpox_fixed.vcf"
REF="ref/NC_063383.fasta"
FAI="ref/NC_063383.fasta.fai"
GFF=$(ls ref/prokka_out/*.gff | head -n 1)
CLEAN_GFF="clean_annotation.gff"

# --- 2. STEP 0: STANDARDIZE FILES ---
echo "--- Step 0: Fixing VCF Header & Standardizing GFF ---"
# Fix VCF header so it matches the reference index
bcftools reheader -f "$FAI" "$VCF" -o "$FIXED_VCF"

# Get the official chromosome name from the .fai file
REAL_NAME=$(head -n 1 "$FAI" | cut -f1)

# Strip FASTA from GFF and fix the chromosome naming mismatch
sed '/^##FASTA/,$d' "$GFF" > temp_clean.gff
PROKKA_NAME=$(grep -v "^#" temp_clean.gff | head -n 1 | cut -f1)

echo "Mapping Prokka ID [$PROKKA_NAME] to Official ID [$REAL_NAME]..."
# Using comma as delimiter in sed to handle pipes in Prokka names
sed "s,${PROKKA_NAME},${REAL_NAME},g" temp_clean.gff > "$CLEAN_GFF"
rm temp_clean.gff

# --- 3. STEP 1-3: VCF TO BED & CONTEXT EXTRACTION ---
echo "--- Step 1-3: Extracting SNP Context (3-bp) ---"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$FIXED_VCF" | \
awk -v name="$REAL_NAME" '{printf "%s\t%d\t%d\t%s/%s\n", name, $2-1, $2, $3, $4}' > MPXV_SNPs.bed

bedtools slop -i MPXV_SNPs.bed -g "$FAI" -b 1 > MPXV_SNPs_3bp.bed
bedtools getfasta -fi "$REF" -bed MPXV_SNPs_3bp.bed -tab -fo MPXV_SNPs_3bp.tab

# --- 4. STEP 5-6: APOBEC3 FILTERING ---
echo "--- Step 5-6: Filtering APOBEC3-like Mutations ---"
paste MPXV_SNPs.bed MPXV_SNPs_3bp.tab | awk '{print $1, $2, $3, $4, $6}' | \
awk '{
    split($4, a, "/"); ref=a[1]; alt=a[2]; context=toupper($5);
    is_apobec=0;
    # Motif: TC->TT (Forward) or GA->AA (Reverse)
    if (context ~ /^TC/ && ref=="C" && alt=="T") is_apobec=1;
    if (context ~ /GA$/ && ref=="G" && alt=="A") is_apobec=1;
    if (is_apobec == 0) print $1"\t"$2"\t"$3;
}' > MPXV_SNPs_noAPOBEC.bed

# --- 5. STEP 7-10: HOTSPOT ANALYSIS ---
echo "--- Step 7-10: Sliding Window Hotspot Analysis (500bp) ---"
bedtools makewindows -g "$FAI" -w 500 > MPXV_windows_500.bed
bedtools coverage -a MPXV_windows_500.bed -b MPXV_SNPs_noAPOBEC.bed > MPXV_window_SNPdensity_500.tab

# Identify Hotspots (Density > 0.01 SNPs/bp)
awk '$4/($3-$2) > 0.01 {print $0 "\tDensity:" $4/($3-$2)}' MPXV_window_SNPdensity_500.tab > MPXV_hotspots_500.tab

# --- 6. STEP 11: ANNOTATION ---
echo "--- Step 11: Annotating Hotspots with Prokka Genes ---"
# Intersect windows with GFF and extract ID/Product attributes
bedtools intersect -a MPXV_hotspots_500.tab -b "$CLEAN_GFF" -wa -wb | \
awk '$0 ~ "\tCDS\t"' | \
perl -ne 'if(/(ID=[^;]+).*?(product=[^;\n]+)/){ 
    @f=split(/\t/); 
    print join("\t", @f[0..2], @f[3], @f[7], $1, $2)."\n"; 
}' | column -t > MPXV_hotspots_annotated.txt

# --- 7. FINAL SUMMARY ---
echo "--- Summary: Top 5 Mutated Gene Products ---"
if [ -s MPXV_hotspots_annotated.txt ]; then
    cut -f 6- MPXV_hotspots_annotated.txt | sort | uniq -c | sort -nr | head -n 5
else
    echo "Check failed: No overlaps found between hotspots and CDS regions."
fi

echo "Pipeline Complete. View results in MPXV_hotspots_annotated.txt"
