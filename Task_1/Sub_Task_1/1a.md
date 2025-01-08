# Subtask 1a - Coverage Analysis

## Objective

In this subtask, we will analyze the coverage of CpG sites across various tissues. We will calculate the coverage as the number of occurrences of each CpG site for each tissue and compute the median and coefficient of variation (CV) for each tissue.



---

## Data Overview

The data used in this analysis is stored in the `PupilBioTest_PMP_revA_task1.csv` file, which contains the following columns:

- `Tissue`: Tissue type (Tissue #1 or Tissue #2).
- `CpG_Coordinates`: The coordinates of the CpG sites.
- `Strand` : Indicates the DNA strand (‘f’ or ‘r’).
- `Methylation Status` : Eight possible patterns (‘000’ to ‘111’).
- `Sample ID` : Unique identifier for each sample.
- `Replicate` : Indicates technical replicates. 



---

## Approach

1. Load the necessary libraries in R.
2. Read the CSV file containing the data.
3. Group the data by tissue and CpG coordinates.
4. Calculate the coverage (the number of occurrences of each CpG site for each tissue).
5. Compute the median and coefficient of variation (CV) for each tissue.



---

## Code Implementation

The following R code demonstrates how to load the data, calculate the coverage, and compute the median and CV for each tissue:

```r
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
```
## Ouput

The output from the above R code will be a table summarizing the median coverage and CV for each tissue, as shown below:
| Tissue   | Median | CV   |
|----------|--------|------|
| Islet    | 64     | 23.2 |
| cfDNA    | 190    | 20.4 |

This output shows that:

The median coverage for Islet is 64, with a CV of 23.2%.
The median coverage for cfDNA is 190, with a CV of 20.4%.

These results suggest that:
The cfDNA tissue has higher coverage with slightly lower variability compared to the Islet tissue.


## Discussion
A high CV indicates a larger variability in coverage across the CpG sites within a tissue.
A low CV suggests that the CpG sites are covered more uniformly across the tissue.
This analysis will be useful for identifying tissues where the sequencing coverage is more consistent, and further experiments can be planned accordingly.

## Conclusion
In this subtask, we successfully analyzed the coverage of CpG sites across tissues and calculated the median and CV for each tissue. This is an essential step in understanding the quality and consistency of the sequencing data.
