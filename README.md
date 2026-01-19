# De Novo *E. coli* Genome Assembly Optimization Study

**Scope:** *E. coli* Genome Assembly Using Velvet & Oases
"Velvet is a de novo genomic assembler specially designed for short-read sequencing technologies. Oases is a de novo transcriptome assembler designed to work on top of Velvet's output, primarily for RNA-Seq data."

## 🧰 Tech Stack & Tools

![Linux](https://img.shields.io/badge/Linux-000?logo=linux\&logoColor=white)
![Bash](https://img.shields.io/badge/Bash-121011?logo=gnu-bash\&logoColor=white)
![Velvet](https://img.shields.io/badge/Velvet-Genome%20Assembler-blue)
![Oases](https://img.shields.io/badge/Oases-Transcriptome%20Assembler-green)
![QUAST](https://img.shields.io/badge/QUAST-Assembly%20QC-orange)
![Python](https://img.shields.io/badge/Python-3776AB?logo=python\&logoColor=white)
![HPC](https://img.shields.io/badge/HPC-Cluster-red)

---

## 2. Objectives

* Assemble the *E. coli* genome using **Velvet** and **Oases**
* Evaluate assemblies at k‑mer sizes **43, 55, 67, 78, and 101**
* Compare performance metrics across k‑mers and assemblers
* Identify the optimal k‑mer size based on assembly statistics

---

## 3. Tech Stack & Tools

![Linux](https://img.shields.io/badge/Linux-000?logo=linux\&logoColor=white)
![Bash](https://img.shields.io/badge/Bash-121011?logo=gnu-bash\&logoColor=white)
![Velvet](https://img.shields.io/badge/Velvet-Genome%20Assembler-blue)
![Oases](https://img.shields.io/badge/Oases-Transcriptome%20Assembler-green)
![QUAST](https://img.shields.io/badge/QUAST-Assembly%20QC-orange)
![Python](https://img.shields.io/badge/Python-3776AB?logo=python\&logoColor=white)
![HPC](https://img.shields.io/badge/HPC-Cluster-red)

---

## 4. Input Data

* `SRR21904868_1.fastq`
* `SRR21904868_2.fastq`

Paired‑end Illumina DNA sequencing reads downloaded from NCBI SRA.

---

## 5. Software Requirements

### Conda Environment Dependencies

* velvet
* oases
* quast
* sra-tools

---

## 6. Workflow Summary

1. Download sequencing data from SRA
2. Set up Conda environment
3. Run Velvet genome assemblies for multiple k‑mers
4. Evaluate Velvet assemblies using QUAST
5. Run Oases assemblies using Velvet outputs
6. Evaluate Oases assemblies using QUAST
7. Compare and interpret results

---

## 7. Step‑by‑Step Methodology

### Step 1: Downloading Sequencing Data

Sequencing reads were downloaded from NCBI SRA using the **SRA Toolkit** on the IU HPC system.

```bash
conda install -c bioconda sra-tools
fasterq-dump SRR21904868
```

This generates paired‑end FASTQ files:

* `SRR21904868_1.fastq`
* `SRR21904868_2.fastq`

---

### Step 2: Conda Environment Setup

A dedicated Conda environment was created to ensure reproducibility.

```bash
conda create --name genome_assembly -c bioconda velvet oases quast
conda activate genome_assembly
```

---

### Step 3: Genome Assembly Using Velvet

Velvet assemblies were performed using k‑mer sizes **43, 55, 67, 78, and 101**.

#### Slurm Script (Velvet)

```bash
#!/bin/bash
#SBATCH --job-name=velvet_assembly
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=msiddhe@iu.edu
#SBATCH --output=velvet_%j.out
#SBATCH -A c01064

source /N/u/msiddhe/miniconda3/bin/activate genome_assembly

KMER_VALUES=(43 55 67 78 101)
READS_DIR=/N/u/msiddhe/Quartz/HPC_Assignment1
OUTPUT_DIR=/N/u/msiddhe/Quartz/HPC_Assignment1/output

mkdir -p $OUTPUT_DIR

for KMER in "${KMER_VALUES[@]}"
do
  mkdir -p $OUTPUT_DIR/velvet_k$KMER
  velveth $OUTPUT_DIR/velvet_k$KMER $KMER -shortPaired -fastq -separate \
    $READS_DIR/SRR21904868_1.fastq \
    $READS_DIR/SRR21904868_2.fastq
  velvetg $OUTPUT_DIR/velvet_k$KMER
done
```

Velvet outputs include:

* `contigs.fa`
* `stats.txt`
* Graph and log files

---

### Step 4: Quality Assessment of Velvet Assemblies

QUAST was used to evaluate all Velvet assemblies together.

```bash
quast velvet_k43/contigs.fa velvet_k55/contigs.fa velvet_k67/contigs.fa \
      velvet_k78/contigs.fa velvet_k101/contigs.fa \
      --min-contig 100 \
      -o quast_results_all_kmers_combined
```

QUAST generates HTML and PDF reports summarizing assembly quality metrics.

---

### Step 5: Assembly Using Oases

Oases was run on Velvet output directories for each k‑mer size.

#### Slurm Script (Oases)

```bash
#!/bin/bash
#SBATCH -A c01064
#SBATCH -J oases_all_kmers
#SBATCH -o oases_all_kmers.out
#SBATCH -e oases_all_kmers.err
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=msiddhe@iu.edu

source /N/u/msiddhe/Quartz/miniconda3/bin/activate genome_assembly

KMER_VALUES=(43 55 67 78 101)
OUTPUT_DIR=/N/u/msiddhe/Quartz/HPC_Assignment1/output

for KMER in "${KMER_VALUES[@]}"
do
  oases $OUTPUT_DIR/oases_k$KMER
done
```

Oases outputs include:

* `transcripts.fa`
* `contigs.fa`
* Graph and log files

---

### Step 6: Quality Assessment of Oases Assemblies

```bash
quast oases_k43/contigs.fa oases_k55/contigs.fa oases_k67/contigs.fa \
      oases_k78/contigs.fa oases_k101/contigs.fa \
      --min-contig 100 \
      -o quast_results_oases_contigs_combined
```

---

### Step 7: Transfer Results to Local Machine

```bash
scp msiddhe@quartz.uits.iu.edu:/path/report.pdf ~/Downloads/
scp msiddhe@quartz.uits.iu.edu:/path/report.html ~/Downloads/
```

---

## 8. Assembly Metrics Explained

* **N50:** Contig length such that 50% of the genome is contained in contigs of this length or longer
* **L50:** Number of contigs contributing to N50
* **GC Content:** Percentage of G/C bases in the assembly

---

## 9. Results and Analysis

### Velvet Assembly Results

| K‑mer | Contigs | Largest Contig (bp) | N50 (bp) |
| ----- | ------- | ------------------- | -------- |
| 43    | 50,840  | 255                 | 128      |
| 55    | 27,558  | 257                 | 120      |
| 67    | 21,086  | 389                 | 147      |
| 78    | 22,388  | 1,248               | 245      |
| 101   | 1,292   | 174,606             | 35,590   |

**Best performance:** k = 101

---

### Oases Assembly Results

Oases produced nearly identical results to Velvet across all k‑mer sizes, including identical optimal performance at **k = 101**.

---

## 10. Comparative Interpretation

* Assembly quality improves consistently with increasing k‑mer size
* Both Velvet and Oases show dramatic improvement from k=78 to k=101
* Oases did **not** significantly alter or improve the Velvet assembly
* GC content remains stable (~50.8%), consistent with *E. coli*

---

## 11. Key Insights

* Optimal k‑mer size for this dataset: **101**
* Final assembly closely matches expected *E. coli* genome size (~4.94 Mb)
* Oases provides no added benefit for genome (DNA‑seq) assembly in this case
* Highlights importance of k‑mer optimization in de novo assembly

---

## 12. Conclusion

This study demonstrates that **k‑mer optimization is critical** for high‑quality de novo genome assembly. For this *E. coli* dataset, a k‑mer size of **101** yielded the best assembly metrics using both Velvet and Oases. Despite Oases being designed for transcriptome assembly, it did not improve genome assembly results beyond Velvet. These findings emphasize the need to match assembly tools appropriately to data type and biological context.
