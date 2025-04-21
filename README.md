## 🧩 Step 1: Gene Identification and Extraction (explained super simply)

### 🎯 **Goal:**
We want to find and save the *exact DNA sequences* of ~80 special genes called **OXPHOS genes** from a big file that contains the entire DNA of a duck (like a giant instruction book).

---

## 🧸 What are we doing and why?

Imagine the duck's genome is a HUGE book with millions of letters (A, T, C, G). Somewhere inside this book are instructions (genes) for making energy. We want to:
1. Make a list of these energy genes (called OXPHOS genes).
2. Use that list to **look up where each gene lives in the big genome book**.
3. **Cut out** just those genes and **save them in a separate file**, so we can study them more easily.

---

## 🧰 What do we need?
- A **list of OXPHOS gene names** (like COX4I1, NDUFA6, ATP5F1A…)
- The **duck genome file** (like a giant FASTA file with all DNA letters)
- The **duck gene map** (a GFF3 or GTF file — tells us where each gene starts and ends)
- Some tools like:
  - `grep`: finds words in files
  - `awk`: smart text cutter
  - `gffread`: extracts gene sequences
  - `bedtools getfasta`: cuts out DNA sequences by coordinates

---

## 🧗 Step-by-Step Instructions

---

### 🟢 **Step 1: Get your list of genes**
You need a list of about 80 genes that are involved in the OXPHOS pathway.

#### 📌 How?
- You can Google: `"OXPHOS genes KEGG pathway"`  
- Or go here: [KEGG Pathway for OXPHOS](https://www.genome.jp/pathway/map00190)
- Write down all the gene names you find into a plain text file called `oxphos_gene_list.txt`.

```
Example content:
NDUFA6
COX4I1
ATP5F1A
...
```
![image](https://github.com/user-attachments/assets/6c87540a-9426-4f05-b914-1420dbcf3844)
--- gene list downloaded and uploaded to oxphos repo as OXPHOS_gene_List.csv
---

### 🟢 **Step 2: Download the duck genome and annotation**

Go to the NCBI genome browser and download:

- **Genome FASTA file** (DNA letters):  
  `GCF_047663525.1_ASM476635v1_genomic.fna.gz`

- **Annotation GFF3 file** (gene map):  
  `GCF_047663525.1_ASM476635v1_genomic.gff.gz`

You can do this using `wget` on the command line:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/663/525/GCF_047663525.1_ASM476635v1/GCF_047663525.1_ASM476635v1_genomic.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/663/525/GCF_047663525.1_ASM476635v1/GCF_047663525.1_ASM476635v1_genomic.gff.gz
```

Then unzip them:

```bash
gunzip GCF_047663525.1_ASM476635v1_genomic.fna.gz
gunzip GCF_047663525.1_ASM476635v1_genomic.gff.gz
```

---

### 🟢 **Step 3: Find the genes in the annotation file**

Let’s say you have a gene called `COX4I1`. You want to find it in the annotation file.

Run:

```bash
grep -i "COX4I1" GCF_047663525.1_ASM476635v1_genomic.gff
```

This will give you lines like:

```
scaffold_1269   RefSeq  gene    19340   20472   .   +   .   ID=gene-COX4I1;Name=COX4I1
```

This tells you the gene lives on `scaffold_1269` from position 19340 to 20472.

---

### 🟢 **Step 4: Extract the DNA sequence of the gene**

There are two easy ways to do this:

#### Option 1: Use `bedtools`

First, make a BED file with all gene coordinates like this:

```
scaffold_1269   19339   20472   COX4I1
```

(Remember: BED files start counting from 0, so we subtract 1 from the start position.)

Save that as `oxphos_genes.bed`.

Now extract the DNA:

```bash
bedtools getfasta -fi GCF_047663525.1_ASM476635v1_genomic.fna -bed oxphos_genes.bed -fo oxphos_genes.fasta
```

You’ll get a FASTA file with your gene sequences!

---

#### Option 2: Use `gffread`

You can extract all genes at once if your GFF file has `gene_name` or `ID=` tags.

```bash
gffread -w oxphos_genes.fasta -g GCF_047663525.1_ASM476635v1_genomic.fna GCF_047663525.1_ASM476635v1_genomic.gff
```

Then you can just open `oxphos_genes.fasta` and grab the sequences for your 80 genes.

---

### 🟢 **Step 5: Double-check your results**

Open the FASTA file in a text editor or use `less`:

```bash
less oxphos_genes.fasta
```
Would you like a ready-made Bash script to automate these steps? Or help with visualizing where these genes are on the genome?
