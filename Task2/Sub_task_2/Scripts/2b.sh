#!/bin/bash

# Exit script on any command failure
set -e

# Define variables
REFERENCE_GENOME="GRCh38.fasta"
TUMOR_R1="PA221MH-lib09-P19-Tumor_S2_L001_R1_001.fastq.gz"
TUMOR_R2="PA221MH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz"
NORM_R1="PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq.gz"
NORM_R2="PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq.gz"

# Step 1: Index the reference genome
echo "Indexing the reference genome..."
bwa index $REFERENCE_GENOME

echo "Reference genome indexed successfully."

# Step 2: Align sequencing reads
echo "Aligning tumor sequencing reads..."
bwa mem $REFERENCE_GENOME $TUMOR_R1 $TUMOR_R2 > Tumor_aligned.sam
echo "Aligning normal sequencing reads..."
bwa mem $REFERENCE_GENOME $NORM_R1 $NORM_R2 > Norm_aligned.sam

echo "Alignment completed."

# Step 3: Convert SAM to BAM and sort
echo "Converting Tumor SAM to sorted BAM..."
samtools view -Sb Tumor_aligned.sam | samtools sort -o Tumor_aligned_sorted.bam
echo "Converting Norm SAM to sorted BAM..."
samtools view -Sb Norm_aligned.sam | samtools sort -o Norm_aligned_sorted.bam

echo "SAM to sorted BAM conversion completed."

# Step 4: Index the BAM files
echo "Indexing Tumor BAM file..."
samtools index Tumor_aligned_sorted.bam
echo "Indexing Norm BAM file..."
samtools index Norm_aligned_sorted.bam

echo "Indexing completed."

# Cleanup intermediate files
echo "Cleaning up intermediate SAM files..."
rm Tumor_aligned.sam Norm_aligned.sam

echo "Pipeline completed successfully! Output BAM and index files are ready."

