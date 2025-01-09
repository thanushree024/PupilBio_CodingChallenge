## Somatic Mutation Identification using Mutect2
This repository provides a workflow for identifying somatic mutations in tumor-normal paired samples using `GATK Mutect2`. Somatic mutations are changes in DNA that occur in cancerous tissue but are absent in normal tissue. This workflow is critical for cancer genomics research and precision oncology.

## Why Mutect2?
Mutect2 is a robust, industry-standard tool for somatic mutation detection. It employs advanced statistical methods to identify mutations with high sensitivity and specificity while effectively handling contamination and sequencing artifacts. This tool is particularly well-suited for:

- Paired Tumor-Normal Analysis: It compares tumor and normal samples to distinguish somatic mutations from germline variants.

- Error Handling: It incorporates features like adaptive error modeling to reduce false positives.

- Customizability: Mutect2 can be easily integrated with downstream annotation tools such as Funcotator for functional analysis of variants.

## Requirements: 
- GATK 4.x or higher
- Reference genome files (`GRCh38.fasta`, `GRCh38.fasta.fai`, `GRCh38.dict`) 
- BAM files (`Tumor_aligned_sorted.bam`, `Norm_aligned_sorted.bam`)

## Mutation Calling Workflow:
This script performs somatic mutation calling and post-processing.

```bash
# Step 1: Set up file paths
TUMOR_BAM="Tumor_aligned_sorted.bam"
NORMAL_BAM="Norm_aligned_sorted.bam"
REFERENCE_GENOME="GRCh38.fasta"
OUTPUT_VCF="somatic_mutations.vcf.gz"

# Step 2: Run Mutect2 for somatic mutation calling
gatk Mutect2 \
    -R $REFERENCE_GENOME \
    -I $TUMOR_BAM \
    -I $NORMAL_BAM \
    --tumor-sample tumor \
    --normal-sample normal \
    -O $OUTPUT_VCF

# Step 3: Filter Mutect2 calls to remove artifacts
FILTERED_VCF="filtered_somatic_mutations.vcf.gz"
gatk FilterMutectCalls \
    -V $OUTPUT_VCF \
    -R $REFERENCE_GENOME \
    -O $FILTERED_VCF

# Step 4: Optional annotation with Funcotator
ANNOTATED_VCF="annotated_somatic_mutations.vcf.gz"
gatk Funcotator \
    -R $REFERENCE_GENOME \
    -V $FILTERED_VCF \
    --output $ANNOTATED_VCF \
    --output-file-format VCF \
    --data-sources-path /path/to/funcotator_dataSources

echo "Somatic mutation calling and filtering complete. Results in: $FILTERED_VCF"
```

### 1. Mutect2: Identifying Somatic Mutations
This step compares the tumor and normal samples to detect somatic mutations while distinguishing them from germline variants.
Mutect2 uses a Bayesian model to evaluate evidence for each variant and classify it as somatic or germline.

### 2. FilterMutectCalls: Removing Artifacts
Raw Mutect2 calls often include sequencing artifacts or low-confidence mutations. This step filters out these false positives, ensuring high-quality results.

### 3. Annotation with Funcotator (Optional)
Funcotator adds functional annotations to variants, providing insights into the biological relevance of detected mutations (e.g., coding vs. non-coding regions, known pathogenic variants).

## Somatic Mutation Identification:
Somatic mutations are mutations that are present in the cancer tissue but absent in the normal tissue. After somatic mutation calling (e.g., using tools like Mutect2, Strelka2, or VarScan2), somatic mutations will appear in the VCF file with specific annotations.

After annotating a VCF file (Variant Call Format file) for somatic mutation identification, you will typically look for certain fields and annotations to distinguish somatic mutations from germline mutations and to estimate background mutations. Here’s how you can interpret the relevant sections in a VCF file:


### Key Columns for Somatic Mutations:

- INFO Field: This field typically contains critical information about the mutation, including its `somatic` status.

Look for annotations like SOMATIC or Somatic in the INFO field, which indicates that the mutation is somatic.
For example, Mutect2 adds a SOMATIC flag to the INFO field for somatic variants.
```javascript
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Mutation is somatic">
```
- FILTER field: The FILTER field may show if a variant passed specific filters or was flagged as somatic. For somatic mutations, the FILTER field could contain values such as `PASS` or `SOMATIC`.
  
- Sample Columns: The individual sample columns (typically named by sample ID) will show the genotype information for each sample. For a somatic mutation:
  - Normal tissue will typically show no variant (or a very low allele frequency, depending on sequencing depth).
  - Cancer tissue will show the presence of the mutation (often with higher allele frequency, depending on the mutation's prevalence in the tumor).

Example of a Somatic Mutation in VCF:
Here is an example of a somatic mutation annotation after running Mutect2:
```CSS
#CHROM  POS     ID        REF  ALT   QUAL  FILTER   INFO
1       123456  .         A    T     100   PASS     DP=100;SOMATIC;VAF=0.35;AF=0.45
```
In this example:
- `SOMATIC`: The SOMATIC tag helps identify those variants that should be considered cancer-specific. In a typical analysis, you compare the normal and cancer tissue VCF files, and variants that are present in the cancer but absent in the normal will be flagged as somatic.
- `VAF`:  A VAF of 0.35 (35%) indicates that 35% of the reads in the cancer sample are supporting the T allele at position 123456, while the rest are supporting the A allele (reference). This suggests that the mutation is present in 35% of the tumor cells, meaning that this mutation is likely not present in all tumor cells, possibly due to tumor heterogeneity or subclonal evolution.
- `AF`:  AF=0.45 suggests that this mutation is found in 45% of the reference population or that it is a common mutation in the population, which can be indicative of a background mutation or germline variant.

## Background Mutation Estimation:
Background mutations are typically the germline mutations (common variants found in both cancer and normal tissues) or the normal genome mutations not associated with cancer.

### To estimate background mutations:

- Germline Variants: These mutations are present in both the normal tissue and the cancer tissue, typically at higher allele frequencies in the normal tissue. These mutations can be identified by comparing the normal tissue VCF file against the cancer tissue VCF file.

- Depth of Coverage (DP): The DP field indicates the sequencing depth at that position. In regions with high background mutations, the depth may be higher than in regions with somatic mutations. This could be a useful indicator of background noise or common polymorphisms.

Example:
```css
#CHROM  POS     ID        REF  ALT   QUAL  FILTER   INFO
1       987654  .         G    A     99    PASS     DP=150;AF=0.05;GERMLINE
```
In this example:
- INFO: `GERMLINE` indicates that this mutation is a part of the individual's normal genetic background and should be considered when estimating the background mutation rate. Since germline mutations are inherited, they should not be included when calculating the somatic mutation burden (mutations specific to cancer). These mutations are part of the natural variation in the population and should be excluded from analyses focused on identifying cancer-specific mutations.
- AF (Allele Frequency): `AF=0.05` suggests that the mutation is common in the population and likely not cancer-specific. Since it is found at a low frequency in the population, it may be a benign variant or normal polymorphism. This helps distinguish it from rare somatic mutations that are typically present in only a subset of cancer cells.
- DP (Depth of Coverage): The `DP=150` indicates that there are 150 sequencing reads at this position, providing evidence for the presence of the G>A mutation. Higher depth of coverage improves the reliability of mutation detection.
