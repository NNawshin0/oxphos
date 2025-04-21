
# 🧬 OXPHOS_Pathway_Altitude_Adaptation

## Overview
This repository documents the pipeline developed to analyze ~80 nuclear-encoded genes in the **oxidative phosphorylation (OXPHOS) pathway** for genetic variation between **high- and low-altitude bird populations**, using the *North American Mallard Genome (NAwild_v1.0)* as reference. The pipeline is inspired by previous work on HIF-pathway genes (e.g., Graham & McCracken 2019, 2021) and adapted to the OXPHOS system.

## 📂 Repository Contents

- `for_loop_trim_filter.txt` – Trimming/filtering script (FASTX Toolkit)
- `for_loop_align_assemble.txt` – BWA/SAMtools-based mapping and sorting
- `remove_orphans_pe.pl` – Perl script to remove unmatched paired-end reads
- `HIF_reference.fasta` – FASTA reference file (example from HIF project)
- `bcf_consensus_test.sh` – BCFtools-based consensus variant script
- `README.md` – This file

## 🧪 Project Workflow

### 🔹 1. Gene Identification and Extraction
**Objective:** Identify and extract full gene models (including UTRs/exons/introns) for ~80 nuclear-encoded OXPHOS genes.

**Steps:**
- Use **KEGG**, **NCBI Gene**, or literature to curate a list of OXPHOS genes.
- Download the *GCF_047663525.1_ASM476635v1* reference genome and annotation files.
- Use `grep`, `awk`, or GFF3 parsers (e.g., `gffread`, BEDTools) to extract coordinates and sequences of each gene.

### 🔹 2. Design of 120mer Hybridization Probes
**Objective:** Design 120 bp oligonucleotide probes to selectively enrich OXPHOS genes for sequencing.

**Steps:**
- Soft-mask extracted gene sequences using **RepeatMasker**.
- Design **120mer probes** at 2× tiling density using in-house scripts or services like Arbor Biosciences (formerly MYcroarray).
- Filter out non-specific probes with genome-wide in silico hybridization tests.

### 🔹 3. Target Enrichment & Sequencing
**Objective:** Capture and sequence OXPHOS gene regions from genomic libraries.

**Steps:**
- Prepare dual-indexed libraries.
- Hybridize with probe set, purify using magnetic beads.
- Sequence using Illumina HiSeq (100 bp paired-end reads, 250–300 bp insert size).

### 🔹 4. Quality Control of Raw Reads
**Objective:** Remove adapters, low-quality reads, and contaminants.

**Script:** `for_loop_trim_filter.txt`

```bash
fastx_clipper -a AGATCGGAAGAGC -l 10 -n -v -i sample_1.fastq -o sample_1_trim.fastq
fastq_quality_filter -v -q 20 -p 90 -i sample_1_trim.fastq -o sample_1_filter.fastq
```

### 🔹 5. Remove Orphaned Read Pairs
**Objective:** Match paired reads and exclude singletons.

**Script:** `remove_orphans_pe.pl`

- Produces: `_sorted.fastq` files and `_singletons.fastq`

### 🔹 6. Sequence Alignment and Sorting
**Objective:** Map quality-filtered reads to reference genes.

**Script:** `for_loop_align_assemble.txt`

```bash
bwa mem HIF_reference.fasta sample_1_sorted.fastq sample_2_sorted.fastq > sample.sam
samtools view -bS sample.sam | samtools sort -o sample.sorted.bam
```

### 🔹 7. Generate VCF from All Individuals
**Objective:** Identify SNPs and small variants across populations.

**Command (example):**

```bash
samtools mpileup -uf OXPHOS_reference.fasta *.sorted.bam | bcftools call -vmO z -o oxphos_all.vcf.gz
```

### 🔹 8. Variant Filtering and Population Structure Analysis
**Objective:** Filter and analyze variants for signs of selection.

**Steps:**
- Decompress and inspect `.vcf.gz` files.
- Use `VCFtools`, `PLINK`, or custom R scripts to calculate:
  - **FST** (between populations)
  - **π (nucleotide diversity)**
  - **Tajima’s D**
- Visualize FST across the genome to identify outlier loci.

### 🔹 9. Optional: Test for Adaptive Introgression
**Objective:** Identify introgressed loci using D-statistics (ABBA–BABA).

**Tool:** [Dsuite](https://github.com/millanek/Dsuite)

```bash
Dsuite Dtrios OXPHOS_tree.nwk oxphos.vcf
```

## 📘 Notes

- All tools should be installed and in `$PATH`. Conda environments recommended.
- All `.vcf` outputs are compressed and must be unzipped or viewed with `zcat`, `bcftools view`, etc.
- Ensure sample naming conventions are consistent across all steps.

## 📚 Citation
If you use this pipeline, please cite:
- Graham & McCracken (2019). *Heredity* 122:819–832
- Graham et al. (2021). *Heredity* 127:107–123
