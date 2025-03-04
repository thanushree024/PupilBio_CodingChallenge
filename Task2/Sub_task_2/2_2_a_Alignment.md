## Alignment and Mutation Calling

Alignment and mutation calling are crucial steps in Next-Generation Sequencing (NGS) analysis. Alignment involves mapping raw sequencing reads to a reference genome to identify their genomic location, typically using tools like BWA or Bowtie. This step ensures that the reads are properly positioned for downstream analysis.

Mutation calling, on the other hand, involves identifying genomic variations such as single nucleotide variants (SNVs), insertions, deletions (indels), or structural variants by comparing aligned reads against the reference genome. Tools like GATK or FreeBayes are commonly used for mutation calling, enabling the detection of disease-associated or functional genomic changes. These steps are foundational in understanding genetic variations in research fields like cancer genomics, population genetics, and precision medicine.

## Data Overview

File Format: *FASTQ*

Samples:
1. `PA221MH-lib09-P19-Tumor_S2_L001_R1_001.fastq.gz`
2. `PA221MH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz`
3. `PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq.gz`
4. `PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq.gz`

Tool Used: `BWA` 

Reference Genome: `GRCh38.fasta`

## Steps to align samples
### 1. Prepare the reference genome:
Indexing the reference genome is crucial because it creates a data structure that allows rapid lookup of sequences during alignment. Without indexing, alignment would be computationally intensive and slow.

#### Command: 
```bash
bwa index GRCh38.fasta
```
### 2. Aligning sequencing reads:
The raw sequencing reads (FASTQ files) are aligned to the reference genome to determine where each read originates.
Paired-end reads (-1 and -2) are used for higher confidence in mapping because they provide positional information from both ends of the fragment.

#### Command: 
```bash
bwa mem GRCh38.fasta PA221MH-lib09-P19-Tumor_S2_L001_R1_001.fastq.gz PA221MH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz > Tumor_aligned.sam
bwa mem GRCh38.fasta PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq.gz PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq.gz > Norm_aligned.sam
```
### 3. Convert SAM to BAM and sort:
The SAM (Sequence Alignment/Map) format is a human-readable output of the alignment, but it is large and inefficient for downstream processing. Converting to BAM (Binary Alignment/Map) format reduces file size and increases processing speed.
Sorting the BAM file arranges the reads by their genomic coordinates, which is essential for tools that perform downstream analyses like variant calling.

#### Command: 
```bash
samtools view -Sb Tumor_aligned.sam | samtools sort -o Tumor_aligned_sorted.bam
samtools view -Sb Norm_aligned.sam | samtools sort -o Norm_aligned_sorted.bam
```

### 4. Index the BAM File
Indexing the sorted BAM file allows rapid retrieval of specific genomic regions during analysis. For example, variant callers need to quickly access reads covering specific positions in the genome.

#### Command: 
```bash
samtools index Tumor_aligned_sorted.bam
samtools index Norm_aligned_sorted.bam
```

The result of this pipeline is a sorted and indexed BAM file (`sample_aligned_sorted.bam`). This file is compact, organized, and ready for downstream tasks like mutation calling, coverage analysis, or visualization in genome browsers.
This workflow ensures efficient processing and prepares the data for accurate variant detection, making it a foundational step in NGS analysis.


