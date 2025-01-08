# Load necessary library
library(dplyr)

# Load the dataset from a CSV file
methylation_data <- read.csv("PupilBioTest_PMP_revA.csv")

# Check the first few rows of the dataset to confirm it's loaded correctly
G#head(methylation_data)

# Check the column names to ensure correct references
colnames(methylation_data)

# Calculate the total coverage for each row (sum of methylation columns)
methylation_data$Single_CpG_Coverage <- rowSums(methylation_data[, c('X.000', 'X.001', 'X.010', 'X.011', 'X.100', 'X.101', 'X.110', 'X.111')])

# Group by tissue type and calculate the median and CV of single CpG coverage
coverage_stats <- methylation_data %>%
  group_by(Tissue) %>%
  summarise(
    Median_CpG_Coverage = median(Single_CpG_Coverage),
    CV_CpG_Coverage = (sd(Single_CpG_Coverage) / mean(Single_CpG_Coverage)) * 100
  )

# Print the results
print(coverage_stats)
