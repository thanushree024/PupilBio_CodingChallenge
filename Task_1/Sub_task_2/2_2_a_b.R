# Load necessary libraries
library(dplyr)
library(tidyr)

# Load the dataset (update the path accordingly)
data <- read.csv("PupilBioTest_PMP_revA.csv")

# Step 1: Perform statistical test for tissue-specific PMPs
# Summarize methylation patterns by Tissue and CpG_Coordinates
methylation_summary <- data %>%
  group_by(Tissue, CpG_Coordinates) %>%
  summarise(across(`X.000`:`X.111`, sum)) %>%
  pivot_longer(cols = `X.000`:`X.111`, names_to = "Pattern", values_to = "Count") %>%
  filter(Count > 0)  # Keep non-zero counts

# Create a contingency table for statistical testing
contingency_tables <- methylation_summary %>%
  group_by(CpG_Coordinates, Pattern) %>%
  summarise(Tissue_1_Count = sum(ifelse(Tissue == "cfDNA", Count, 0)),
            Tissue_2_Count = sum(ifelse(Tissue == "Islet", Count, 0))) %>%
  ungroup() %>%
  mutate(Total_Count = Tissue_1_Count + Tissue_2_Count)

# Apply Chi-Square or Fisher's Exact Test for each PMP
results <- contingency_tables %>%
  rowwise() %>%
  mutate(
    p_value = ifelse(Total_Count < 5,
                     fisher.test(matrix(c(Tissue_1_Count, Tissue_2_Count, 
                                          Total_Count - Tissue_1_Count, Total_Count - Tissue_2_Count),
                                        nrow = 2))$p.value,
                     chisq.test(matrix(c(Tissue_1_Count, Tissue_2_Count, 
                                         Total_Count - Tissue_1_Count, Total_Count - Tissue_2_Count),
                                       nrow = 2))$p.value)
  ) %>%
  arrange(p_value)

# Filter for PMPs with significant differentiation (e.g., p-value < 0.05)
significant_pmps <- results %>%
  filter(p_value < 0.05)

# Step 2: Calculate Mean Variant Read Fraction (VRF)
vrf <- data %>%
  group_by(Tissue, CpG_Coordinates) %>%
  summarise(across(`X.000`:`X.111`, mean)) %>%
  pivot_longer(cols = `X.000`:`X.111`, names_to = "Pattern", values_to = "Mean_VRF")

# Merge VRF results with significant PMPs
final_results <- significant_pmps %>%
  left_join(vrf, by = c("CpG_Coordinates", "Pattern"))

# Save results to a file
write.csv(final_results, "significant_pmps_with_vrf.csv", row.names = FALSE)

# Output
print("Significant PMPs with VRF saved to 'significant_pmps_with_vrf.csv'")

