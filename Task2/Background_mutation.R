# Load necessary libraries
library(dplyr)

# Load VCF data
# Replace 'normal_sample.vcf' with the actual VCF file path
vcf_data <- read.table("results/somatic_variants_mutect2.vcf.gz", header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)

# Assign column names to the VCF data
colnames(vcf_data) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Normal")

# Function to extract a specific field (e.g., AF or DP) from the INFO column
extract_field <- function(info, field) {
  matches <- regmatches(info, regexpr(paste0(field, "=\\d+(\\.\\d+)?"), info))
  if (length(matches) > 0) {
    return(as.numeric(sub(paste0(field, "="), "", matches)))
  } else {
    return(NA)
  }
}

# Extract allele frequency (AF) values from the INFO column
vcf_data$AF <- sapply(vcf_data$INFO, extract_field, field = "AF")

# Calculate the median background mutation level
median_background_mutation <- median(vcf_data$AF, na.rm = TRUE)
cat(sprintf("Median Background Mutation Level: %.6f\n", median_background_mutation))

# Define confidence threshold (e.g., 5% or 0.05)
confidence_threshold <- 0.05

# Calculate the minimum reads required for confident mutation calling
min_reads_required <- confidence_threshold / median_background_mutation
cat(sprintf("Minimum Reads Required: %.2f\n", min_reads_required))

head(vcf_data$AF)
summary(vcf_data$AF)
