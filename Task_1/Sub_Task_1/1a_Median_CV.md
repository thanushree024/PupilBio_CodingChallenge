# Subtask 1a - Coverage Analysis

## Objective

In this subtask, we will perform a coverage analysis on a dataset containing phased methylation patterns (PMPs) across different tissue types. Specifically, it calculates the single CpG coverage for each sample and computes the median and coefficient of variation (CV) for the coverage of CpG sites in various tissue types.



---

## Data Overview

The data used in this analysis is stored in the `PupilBioTest_PMP_revA_task1.csv` file, which contains the following columns:

- `Tissue`: Tissue type (Tissue #1 or Tissue #2).
- `CpG_Coordinates`: The coordinates of the CpG sites.
- `Strand` : Indicates the DNA strand (‘f’ or ‘r’).
- `Methylation Status` : Eight possible patterns (‘000’ to ‘111’).
- `Sample ID` : Unique identifier for each sample.
- `Replicate` : Indicates technical replicates. 

## Approach
1. Load the necessary libraries in R
2. Loading the dataset
3. Inspecting the dataset
4. Checking column names
5. Calculating the single CpG coverage
6. Grouping Data by Tissue Type and Calculating Statistics
7. Output


## Code Implementation

This R code demonstrates how to load the data, calculate the Single CpG coverage, and compute the median and CV for each tissue.

```r
# Load necessary library
library(dplyr)

# Load the dataset from a CSV file
methylation_data <- read.csv("PupilBioTest_PMP_revA.csv")

# Check the first few rows of the dataset to confirm it's loaded correctly
#head(methylation_data)

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

```

## Ouput

The output from the above R code will be a table summarizing the median coverage and CV for each tissue, as shown below:
| Tissue   | Median_CpG_Coverage | CV_CpG_Coverage   |
|----------|--------|------|
| Islet    | 84    | 114 |
| cfDNA    | 484   | 132 |
 
## Explanation of Median and CV in Coverage Analysis
### Median
Median is the middle value in a set of data when the values are sorted in order. In coverage analysis, it is used to summarize the central tendency of the coverage data for each tissue. The median provides a robust measure of central location, which is less sensitive to outliers and skewed distributions compared to the mean. For example, if one CpG site has significantly higher coverage, the median still gives a good representation of the typical coverage across all sites.

### Coefficient of Variation (CV)
CV is the ratio of the standard deviation to the mean, often expressed as a percentage. It is used to measure the relative variability in coverage between tissues. A higher CV indicates more variability, while a lower CV suggests more consistent coverage across the CpG sites within a tissue. It is particularly useful when comparing the variability of coverage across tissues with different median values.

## Conclusion
*Median* methylation counts differ between cfDNA and Islet tissues, highlighting distinct epigenetic profiles and supporting the hypothesis that PMPs provide specificity for tissue differentiation.

*Coefficient of variation (CV)* demonstrates that cfDNA has higher variability in methylation patterns compared to Islet tissue. This may reflect biological heterogeneity in cfDNA due to its multi-tissue origin, whereas Islet tissue shows more consistent methylation patterns.
These findings emphasize the potential of PMPs as biomarkers but also highlight the importance of addressing variability to ensure reliability.
