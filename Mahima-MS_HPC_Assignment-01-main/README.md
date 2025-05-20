# Mahima-MS_Precision Medicine_Assignment-01

Scope: E. coli Genome Assembly Using Velvet & Oases
"Velvet is a de novo genomic assembler specially designed for short read sequencing technologies. Oases is a de novo transcriptome assembler designed to work on top of Velvet's output, primarily for RNA-Seq data."

Assignment Description: This assignment aims to assemble the Escherichia coli genome using short-read Illumina DNA sequencing data (SRA identifier: SRR21904868). Two de-novo genome assembly tools, Velvet and Oases, were used to perform the assembly with varying k-mer sizes. The different k-mers  were taken, and a comparison of the results between Velvet and Oases was conducted, highlighting key differences and performance insights.

Objective:
-> Assemble the E. coli genome using Velvet and Oases.
-> Results for different k-mer sizes (43, 55, 67, 78, 101).
-> Comparing the results and performance of Velvet and Oases.
-> Providing a detailed step-by-step analysis.

Programming Language:
Bash
Python (for running QUAST) ---> Miniconda

Input Files: 
SRR21904868_1.fastq
SRR21904868_2.fastq


Required Files/Dependencies:
Oases installed in a conda environment.
Velvet installed in the same environment for preprocessing.
Short-read FASTQ files.


Required Libraries:
Velvet
Oases
QUAST


# Execution Steps:

# Step 1: Downloading the Escherichia coli Genome
We intialize the first process by downloading short-read illumina DNA sequencing data from the SRA database using the SRA toolkit. We initiate this process by connecting it to the IU HPC System through Quartz server. 

Install SRA Toolkit using conda:
bashscript:
conda install -c bioconda sra-tools

Download the Genome Data: To download the sequencing data SRR21904868 in FASTQ format, use following command:
bashscript:
fasterq-dump SRR21904868
This will download paired-end reads for Escherichia coli as SRR21904868_1.fastq and SRR21904868_2.fastq.


# Step 2: Set Up Miniconda Environment
It is necessary to set up miniconda environment so that it can include all necessary tools required for the analysis.

Create Conda Environment:
bashscript:
conda create --name genome_assembly -c bioconda velvet oases quast

Activate the Environment:
bashscript:
conda activate genome_assembly


# Step 3: Analysis using Velvet (De novo genome assembly)
k-mers sizes: 43, 55, 67, 78, 101
Run Velvet with different k-mer sizes as mentioned. The step-by-step process is as follows:


a. Prepare Slurm Script:
Create a Slurm script for all the different k-mers size for job submission.

bashscript:
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

Activate Conda environment
source /N/u/msiddhe/miniconda3/bin/activate genome_assembly

Set k-mer values and paths
KMER_VALUES=("43" "55" "67" "78" "101")
READS_DIR=/N/u/msiddhe/Quartz/HPC_Assignment1
OUTPUT_DIR=/N/u/msiddhe/Quartz/HPC_Assignment1/output

mkdir -p $OUTPUT_DIR

for KMER in "${KMER_VALUES[@]}"
do
    mkdir -p $OUTPUT_DIR/velvet_k$KMER
    velveth $OUTPUT_DIR/velvet_k$KMER $KMER -shortPaired -fastq -separate $READS_DIR/SRR21904868_1.fastq $READS_DIR/SRR21904868_2.fastq
    velvetg $OUTPUT_DIR/velvet_k$KMER
done

echo "All k-mer processing completed!"


b. Submit Job to Slurm:
bashscript:
sbatch velvet_all_kmers.slurm

Velvet Output Files: After Velvet runs successfully, we will receive the mail that the run is completed and followed by it will produce an output directory (e.g., velvet_k43, velvet_k55) containing files like: contigs.fa: The FASTA file of contigs, for all the k-mers.
Velvet Output (FASTA files for contigs):
contigs.fa
Graph
LastGraph
Log
PreGraph
Roadmaps
Sequences
stats.txt

Run QUAST on Velvet Output:
bashscript:
quast /N/u/msiddhe/Quartz/HPC_Assignment1/output/velvet_k43/contigs.fa \
      /N/u/msiddhe/Quartz/HPC_Assignment1/output/velvet_k55/contigs.fa \
      /N/u/msiddhe/Quartz/HPC_Assignment1/output/velvet_k67/contigs.fa \
      /N/u/msiddhe/Quartz/HPC_Assignment1/output/velvet_k78/contigs.fa \
      /N/u/msiddhe/Quartz/HPC_Assignment1/output/velvet_k101/contigs.fa \
      --min-contig 100 \
      -o /N/u/msiddhe/Quartz/HPC_Assignment1/quast_results_all_kmers_combined

Velvet Quast Output Files:
basic_stats/
icarus_viewers/
report.txt
transposed_report.tex
transposed_report.tsv
transposed_report.txt
quast.log
report.html
report.pdf
report.tex
report.tsv


After we get the Quast result we will analyze the report.html/report.pdf file for the different k-mers value and analyze it.
The HTML report shows a heatmap that visually represents the quality assessment of genome assemblies for different k-mer sizes using Velvet.


# Step 4: Oases Assembly (Transcriptome assembly)
Once Velvet has generated the intermediate files, we will run Oases for transcript assembly.

Run Oases with Different k-mer Sizes:
k-mers sizes: 43, 55, 67, 78, 101
Run oases with different k-mer sizes as mentioned. The step-by-step process is as follows:

a. Prepare Slurm Script for Oases:
#!/bin/bash
#SBATCH -A c01064                    # Allocation ID
#SBATCH -J oases_all_kmers            # Job name
#SBATCH -o oases_all_kmers.out        # Standard output
#SBATCH -e oases_all_kmers.err        # Standard error
#SBATCH -t 4:00:00                    # Time limit
#SBATCH -N 1                          # Number of nodes
#SBATCH -n 4                          # Number of cores
#SBATCH --mail-type=BEGIN,END,FAIL    # Notifications
#SBATCH --mail-user=msiddhe@iu.edu    # Your email for notifications

Activate conda environment
source /N/u/msiddhe/Quartz/miniconda3/bin/activate genome_assembly

Array of k-mer values
KMER_VALUES=("43" "55" "67" "78" "101")

Base directories
OUTPUT_DIR=/N/u/msiddhe/Quartz/HPC_Assignment1/output

Run Oases for each k-mer
for KMER in "${KMER_VALUES[@]}"
do
    echo "Running Oases for k-mer $KMER..."
    oases $OUTPUT_DIR/oases_k$KMER
done

echo "Oases assembly for all k-mers completed!"


b. Submit Job to Slurm:
sbatch oases_all_kmers.slurm

Oases Output Files: After Oases completes, we will find files in the respective k-mers directories such as:
contig-ordering.txt
contigs.fa
Graph
Graph2
LastGraph
Log
PreGraph
Roadmaps
Sequences
stats.txt
transcripts.fa


# Step 5: QUAST for Quality Assessment
QUAST (Quality Assessment Tool for Genome Assemblies) will be used to evaluate the assemblies for both Velvet and Oases.
quast /N/u/msiddhe/Quartz/HPC_Assignment1/output/oases_k43/contigs.fa \
      /N/u/msiddhe/Quartz/HPC_Assignment1/output/oases_k55/contigs.fa \
      /N/u/msiddhe/Quartz/HPC_Assignment1/output/oases_k67/contigs.fa \
      /N/u/msiddhe/Quartz/HPC_Assignment1/output/oases_k78/contigs.fa \
      /N/u/msiddhe/Quartz/HPC_Assignment1/output/oases_k101/contigs.fa \
      --min-contig 100 \
      -o /N/u/msiddhe/Quartz/HPC_Assignment1/quast_results_oases_contigs_combined

echo "QUAST analysis for Oases contigs completed!"


Quast Oases Output File:
basic_stats
icarus_viewers
icarus.html
quast.log
report.html
report.pdf
report.txt
report.tex
report.tsv
transposed_report.tex
transposed_report.tsv
transposed_report.txt

After we get the Quast result we will analyze the report.html/report.pdf file for the different k-mers value and analyze it. 
The HTML report shows a heatmap that visually represents the quality assessment of genome assemblies for different k-mer sizes using Oases.


# Step 6: Transfering the report to the local machine:
scp msiddhe@quartz.uits.iu.edu:/N/u/msiddhe/Quartz/HPC_Assignment1/quast_results_all_kmers_combined/report.pdf C:/Users/mmsid/Downloads/
scp msiddhe@quartz.uits.iu.edu:/N/u/msiddhe/Quartz/HPC_Assignment1/quast_results_all_kmers_combined/report.html C:/Users/mmsid/Downloads/
scp msiddhe@quartz.uits.iu.edu:/N/u/msiddhe/Quartz/HPC_Assignment1/quast_results_oases_all_kmers/report.pdf C:/Users/mmsid/Downloads/
scp msiddhe@quartz.uits.iu.edu:/N/u/msiddhe/Quartz/HPC_Assignment1/quast_results_oases_all_kmers/report.html C:/Users/mmsid/Downloads/


--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Analyzing the Velvet and the Oases Report

Note: 
N50: The length of the shortest contig in the set that contains the fewest contigs whose combined length represents at least 50% of the assembly.
L50: The number of contigs equal to or longer than N50.
GC content: The percentage of nitrogenous bases that are either guanine or cytosine.

## Analysis of Velvet Assembly Report:

K-mer size comparison:
K=43: 50,840 contigs, largest 255 bp, N50 128 bp
K=55: 27,558 contigs, largest 257 bp, N50 120 bp
K=67: 21,086 contigs, largest 389 bp, N50 147 bp
K=78: 22,388 contigs, largest 1,248 bp, N50 245 bp
K=101: 1,292 contigs, largest 174,606 bp, N50 35,590 bp

Optimal parameters:
K=101 produces the best results with fewest contigs, largest contig, and highest N50.
Total assembly length (4,940,551 bp) is closest to the expected E. coli genome size.

Performance analysis:
As k-mer size increases, the number of contigs generally decreases and contig size increases.
Significant improvement from k=78 to k=101 in all metrics.
K=101 is the only size producing contigs ≥1000 bp in significant numbers (251).

## Analysis of Oases Assembly Report:

K-mer size comparison:
K=43: 50,840 contigs, largest 255 bp, N50 128 bp
K=55: 27,582 contigs, largest 257 bp, N50 120 bp
K=67: 21,094 contigs, largest 389 bp, N50 147 bp
K=78: 22,387 contigs, largest 1,248 bp, N50 245 bp
K=101: 1,292 contigs, largest 174,606 bp, N50 35,590 bp

Optimal parameters:
K=101 produces the best results, mirroring the Velvet assembly.
Total assembly length (4,940,551 bp) matches the Velvet assembly.

Performance analysis:
The trend of improvement with increasing k-mer size is consistent with Velvet.
K=101 shows dramatic improvement in all metrics compared to smaller k-mer sizes.

## Comparison and Differences:

Similarity:
The results for Velvet and Oases are nearly identical, especially for k=101.
Both assemblers show the same trends in improvement as k-mer size increases.

Minor differences:
Slight variations in contig numbers for k=55 (Velvet: 27,558, Oases: 27,582)
Minimal differences in total length for smaller k-mer sizes
(Note: Difference are highlighted inthe report that is generated for both Velvet and Oases)

Interpretation:
The similarity suggests that Oases, which builds on Velvet's output, did not significantly modify the assembly for this dataset.
This could be due to the nature of the E. coli genome or the specific characteristics of this sequencing data.

## Highlights:
Both Velvet and Oases perform best with k=101 for this E. coli genome assembly.
The optimal assembly (k=101) produces 1,292 contigs with a largest contig of 174,606 bp and N50 of 35,590 bp.
The total assembly length of 4,940,551 bp closely matches the expected E. coli genome size.
There is a significant improvement in assembly quality as k-mer size increases, with a dramatic jump from k=78 to k=101.

Surprisingly, Velvet and Oases produce nearly identical results, suggesting that for this particular dataset, Oases did not substantially improve upon the Velvet assembly.
The GC content remains consistent across different k-mer sizes and between assemblers, around 50.8%, which is typical for E. coli.
This analysis demonstrates the importance of optimizing k-mer size in de novo genome assembly and highlights the similar performance of Velvet and Oases for this specific E. coli dataset.

Oases is typically used for transcriptome assembly and often produces different results from Velvet. Some points to consider:
The similarity in results suggests that for this particular E. coli dataset, Oases didn't significantly modify the Velvet assembly.
This could be due to the nature of the E. coli genome or the specific characteristics of this sequencing data.
It might be worth investigating why Oases didn't produce different results, as it's designed for transcriptome assembly rather than genome assembly. 


## Conclusion:
For this E. coli genome assembly, a k-mer size of 101 produced the best results for both Velvet and Oases. The optimal assembly yielded 1,292 contigs with a largest contig of 174,606 bp and an N50 of 35,590 bp. The total assembly length of 4,940,551 bp closely matches the expected E. coli genome size.
The k-mer size of 101 yielded the best results because it balances specificity and coverage effectively for this E. coli dataset. Larger k-mers improve specificity by minimizing incorrect links between reads, which is especially useful for repetitive genome regions. However, they also demand higher coverage to form connections. In this case, k=101 proved to be the optimal point where the k-mer was large enough to resolve repetitive regions while still ensuring enough coverage for accurate assembly.

For this particular E. coli genome assembly, there doesn't seem to be any major advantage in using Oases over Velvet. The fact that both tools produced nearly identical results indicates that Oases, despite its design for transcriptome assembly, did not enhance the Velvet assembly for this specific genomic dataset. This unexpected result underscores the value of trying different assembly tools, as their effectiveness can vary depending on the nature of the data being processed.

In conclusion for this E. coli dataset, Oases did not provide any significant improvement, which is consistent with the data.



