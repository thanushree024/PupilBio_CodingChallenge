########### Median Bar Plot ##########################################
median_plot <- ggplot(coverage_stats, aes(x = Tissue, y = Median, fill = Tissue)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(
    title = "Median CpG Coverage Across Tissues",
    x = "Tissue",
    y = "Median Coverage"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set2")

############################# CV Bar Plot ######################################

# Bar plot for Coefficient of Variation (CV) across Tissues
cv_plot <- ggplot(coverage_stats, aes(x = Tissue, y = CV, fill = Tissue)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(
    title = "Coefficient of Variation (CV) of CpG Coverage",
    x = "Tissue",
    y = "CV (%)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set2")

# Display plots
median_plot
cv_plot


coverage_stats <- coverage_data %>%
  group_by(Tissue) %>%
  summarise(
    Median = median(Coverage),
    CV = (sd(Coverage) / mean(Coverage)) * 100,
    SD = sd(Coverage)
  )

######################### scatter plot #####################

# Scatter plot with error bars
scatter_plot <- ggplot(coverage_stats, aes(x = Tissue, y = Median)) +
  geom_point(size = 4, color = "blue") +
  geom_errorbar(aes(ymin = Median - SD, ymax = Median + SD), width = 0.2) +
  labs(
    title = "Median CpG Coverage with Variability",
    x = "Tissue",
    y = "Median Coverage"
  ) +
  theme_minimal()

scatter_plot

################################ Heat Map ######################################

library(tidyr)            # Load the package

# Prepare data for Heatmap
heatmap_data <- coverage_stats %>%
  pivot_longer(cols = c(Median, CV), names_to = "Metric", values_to = "Value")

# Heatmap
heatmap_plot <- ggplot(heatmap_data, aes(x = Tissue, y = Metric, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title = "Heatmap: CpG Coverage Statistics Across Tissues",
    x = "Tissue",
    y = "Metric"
  ) +
  theme_minimal()

heatmap_plot

