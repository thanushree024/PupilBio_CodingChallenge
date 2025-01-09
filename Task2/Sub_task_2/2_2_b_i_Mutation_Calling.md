## Somatic Mutation Identification using Mutect2
This repository provides a workflow for identifying somatic mutations in tumor-normal paired samples using `GATK Mutect2`. Somatic mutations are changes in DNA that occur in cancerous tissue but are absent in normal tissue. This workflow is critical for cancer genomics research and precision oncology.

## Why Mutect2?
Mutect2 is a robust, industry-standard tool for somatic mutation detection. It employs advanced statistical methods to identify mutations with high sensitivity and specificity while effectively handling contamination and sequencing artifacts. This tool is particularly well-suited for:

- Paired Tumor-Normal Analysis: It compares tumor and normal samples to distinguish somatic mutations from germline variants.

- Error Handling: It incorporates features like adaptive error modeling to reduce false positives.

- Customizability: Mutect2 can be easily integrated with downstream annotation tools such as Funcotator for functional analysis of variants.

## Requirements: 
- GATK 4.x or higher
- Reference genome files (`.fasta`, `.fasta.fai`, `.dict`)
- BAM files (`tumor.bam`, `normal.bam`, `.bai indices`)

## Data Overview:
