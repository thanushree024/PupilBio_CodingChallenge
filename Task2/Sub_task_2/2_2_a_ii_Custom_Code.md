## Mutation Analysis with VCF Files
This script analyzes genetic variation using VCF (Variant Call Format) files. It counts the number of Single Nucleotide Polymorphisms (SNPs) and Insertions/Deletions (INDELs), calculates the total number of mutations, and estimates the background mutation rate based on the callable genome size. The analysis is done using the vcfR package in R.

Data Overview:
`SNP.vcf`and `Indel.vcf`[Link to VCF File]()

### Installation
If you don't have the vcfR package installed, you can do so by running the following command in R:
```r
install.packages("vcfR")
```
### Code Implementation
```r
# Load necessary library
library(vcfR)

# Load VCF files
snv_file <- "SNVs.vcf"
indel_file <- "Indels.vcf"

# Function to count mutations
count_mutations <- function(vcf_file) {
  vcf_data <- read.vcfR(vcf_file)
  mutation_count <- nrow(vcf_data@fix)  # Count rows in the 'fix' slot, which contains variant information
  return(mutation_count)
}

# Count SNPs and INDELs
snv_count <- count_mutations(snv_file)
indel_count <- count_mutations(indel_file)

# Define callable genome size (e.g., hg38 ~3 billion bases)
callable_genome_size <- 3e9  

# Calculate mutation rate
total_mutations <- snv_count + indel_count
mutation_rate <- total_mutations / callable_genome_size

# Output results
cat(sprintf("Total SNPs: %d\n", snv_count))
cat(sprintf("Total INDELs: %d\n", indel_count))
cat(sprintf("Total Mutations: %d\n", total_mutations))
cat(sprintf("Background Mutation Rate: %.2e mutations per base pair\n", mutation_rate))
```
Input Files: The script uses two VCF files:

`SNVs.vcf` containing SNP data

`Indels.vcf` containing INDEL data

Mutation Counting: The count_mutations function reads the VCF file and counts the number of variants by looking at the fix slot, which contains the variant information.

Background Mutation Rate: The script calculates the background mutation rate by dividing the total number of mutations (SNPs + INDELs) by the callable genome size (typically ~3 billion bases for the human genome).

### Results:
```mathematics
Total SNPs: 256
Total INDELs: 0
Total Mutations: 256
Background Mutation Rate: 8.53e-08 mutations per base pair
```
