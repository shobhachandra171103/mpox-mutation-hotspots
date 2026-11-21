# Mpox Mutation Hotspots & Selection Pressure (`mpox-mutation-hotspots`)

> **Objective 1** of the B.Sc. (Life Sciences) 4th-Year Dissertation: *"Multi-Scale Genomic Architecture of Monkeypox Virus Human Adaptation: Mutation Hotspots, Structural Variants, and Ortholog Divergence in Sustained-Transmission Lineages"* (Department of Zoology, Deshbandhu College, University of Delhi, 2026).

## Overview

This repository documents the computational workflow and analytical outputs for mapping genome-wide single nucleotide polymorphism (SNP) hotspots—including analyses independent of APOBEC3/APOBEC3F signatures—alongside 3-bp trinucleotide context profiling and selection pressure (`dN/dS`) estimation across **1,236 high-quality Monkeypox virus (MPXV) genomes**.

## Dataset & Filtering Criteria

* **Source:** NCBI Virus (`https://www.ncbi.nlm.nih.gov/labs/virus/`), retrieved as of **02 January 2026**.
* **Cohort Size:** **1,236** complete MPXV genomes.
* **Quality Filters:**
  * Minimum sequence length: **170 kb**
  * Nucleotide completeness: **Complete**
  * Ambiguous characters: **None**
* **Reference Genome:** `NC_063383` (plus custom `.gbk` reference databases).
* **Metadata:** `.csv` compilation including Accession ID, collection date, country, geographic location, and manually added latitude/longitude coordinates.

## Methodology Pipeline

1. **FASTA Preparation & Annotation:**
   * Generate individual FASTA files per accession using Python.
   * Annotate all 1,236 genomes using **PROKKA** to extract `.faa` (predicted proteins), `.gbk` (annotated genomes), and `.gff` coordinate files.
   * Combine `.faa` files and collapse redundancy with **CD-HIT** to build a non-redundant (`nr`) protein set.
2. **Phylogeny & SNP Calling:**
   * Determine optimal k-mer size (`k = 19`) with `kchooser` and call SNPs (`.vcf`) / SNP-based Maximum Likelihood (ML) phylogeny using **kSNP4** (with Jellyfish k-mer partitioning), visualized in **Microreact** with geographic color-coding.
   * Align sequences via **MAFFT** and call variants using **SNP-sites** (Galaxy).
3. **3-bp Window Context & APOBEC3 Filtering:**
   * Convert `.vcf` to `.bed` via `bcftools` and `bedtools`.
   * Index reference sequences with `samtools` (`.fai`) to expand SNPs into **3-bp windows** and extract trinucleotide contexts.
   * Filter APOBEC3-compatible motifs (`TC -> TT` forward / `GA -> AA` reverse) to isolate non-APOBEC3 hotspots.
4. **Sliding Window Hotspot Analysis & Gene Mapping:**
   * Run **500 bp sliding window** analysis (`395` windows across the ~197 kb genome) using `bedtools`.
   * Classify windows with **SNP density > 0.01 SNPs/bp** as mutation hotspots (`257` hotspot windows identified).
   * Intersect hotspot coordinates with PROKKA `.gff` annotations using `bedtools intersect`.
5. **Selection Pressure (`dN/dS`) & Multi-Scale Visualization:**
   * Calculate synonymous (`dS`) and non-synonymous (`dN`) substitution ratios (`dN/dS`) in Python (`pandas`).
   * Generate SNP density plots, 3-bp trinucleotide motif charts, top mutated gene rankings, annotated Manhattan plots, R-based semi-log `dN/dS` selection plots, and `pycirclize` Circos plots.

## Key Findings Summary

* **Terminal Hotspot Concentration:** Out of **395** total 500 bp windows, **257 windows** (`density > 0.01`) were identified as hotspots, with peak density reaching **0.738** (`369` SNPs in `190,000–190,500 bp`) concentrated in the Inverted Terminal Repeats (ITRs) and host-range / immune-evasion accessory regions.
* **3-bp Context & Codon Bias:** Trinucleotide motif distributions reflect strong AT-richness (~67% AT) and mutational/selection pressures, with high representation of AT-rich motifs (such as `TCA` and `ACT`) and low counts for GC-heavy motifs (`CCC`, `GGG`).
* **Top Mutated Genes & Selection:** High-mutation loci include hypothetical proteins (`PHDGHHDA_00204`, `PHDGHHDA_00199`, `PHDGHHDA_00208`, `PHDGHHDA_00209`) as well as immune-modulating and host-interaction genes (e.g., ankyrin repeat proteins `OPG023`/`OPG025`, pro-viral epidermal growth factor, `Protein C10`, and soluble interferon alpha/beta receptor `OPG204`). Selection analysis (`114` of `152` genes showing `dN/dS > 1` in the selection dataset) highlights strong positive selection in genes such as `PHDGHHDA_00191` (pseudokinase `OPG198`, ratio `18.0`) and `PHDGHHDA_00118` (major core protein `OPG129`, ratio `11.0`).

## Repository Structure

```text
mpox-mutation-hotspots/
├── README.md
├── scripts/          # Bash, Python, and R scripts for preprocessing, filtering, hotspots, dN/dS, and plotting
├── results/          # Summary tables and generated figures (SNP density, Manhattan, dN/dS, Circos)
└── docs/             # Methodology notes and workflow diagrams
