# Load necessary libraries
library(dplyr)
library(ggplot2)
# Load data
data <- read.csv("PupilBioTest_PMP_revA_task1.csv")
# Calculate coverage (number of occurrences of each CpG site for each tissue)
coverage_data <- data %>%
  group_by(Tissue, CpG_Coordinates) %>%
  summarise(Coverage = n())

# Calculate median and CV for each tissue
coverage_stats <- coverage_data %>%
  group_by(Tissue) %>%
  summarise(
    Median = median(Coverage),
    CV = (sd(Coverage) / mean(Coverage)) * 100
  )

# Print the results
print(coverage_stats)
