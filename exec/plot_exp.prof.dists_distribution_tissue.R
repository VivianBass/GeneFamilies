require(GeneFamilies)
options(mc.cores = getMcCores())
library(dotenv)

library(ggplot2)
library(ggsignif)
library(gridExtra)

library(dplyr)
library(purrr)

library(tidyverse)
library(rstatix)
library(ggpubr)
library(tibble)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

load(file.path(output_data_dir, "expression_profile_distances_statistics.RData"))
load(file.path(output_data_dir, "ExpressionProfileDistances.RData"))

load("data/data.RData")

load(file.path(output_data_dir, "/gene_family_expression_dists.RData"))
load(file.path(output_data_dir, "/boxplot_expression_distances.RData"))

significance_level <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")  
}


load(file.path(output_data_dir, "expression_distances_tissue.RData"))

# Plotting process

# Creating a Boxplot for Mean Distances
# Calculate the appropriate height for the PDF based on the number of facets
n_facets <- length(unique(mean_p.df_long$variables))
plot_height <- 8 * n_facets  # 8 inches per facet

boxplot_mean <- ggplot(mean_p.df_long, aes(x = Type, y = value, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Mean Expression Distances Tissue", y = "Distance") +
  theme_pubr(border = TRUE) +
  geom_signif(comparisons = list(c("Ortholog", "Paralog")), map_signif_level = TRUE) +
  theme(plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, r = 0 ,b = 20, l = 0), hjust = 0.50),
        axis.title.x = element_text(size = 10, margin = margin(t = 20, r = 0, b = 20, l = 0), hjust = 0.50),
        axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20))) +
  facet_grid(variables ~ ., scales = "free_y", space = "free_y")

ggsave(file.path(results_dir, "tissues_mean_boxplot12.pdf"), plot = boxplot_mean, width = 12, height = plot_height)


# Creating a Boxplot for Median Distances
# Calculate the appropriate height for the PDF based on the number of facets
n_facets <- length(unique(median_p.df_long$variables))
plot_height <- 8 * n_facets  # 8 inches per facet

boxplot_mean <- ggplot(median_p.df_long, aes(x = Type, y = value, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Median Expression Distances Tissue", y = "Distance") +
  theme_pubr(border = TRUE) +
  geom_signif(comparisons = list(c("Ortholog", "Paralog")), map_signif_level = TRUE) +
  theme(plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, r = 0 ,b = 20, l = 0), hjust = 0.50),
        axis.title.x = element_text(size = 10, margin = margin(t = 20, r = 0, b = 20, l = 0), hjust = 0.50),
        axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20))) +
  facet_grid(variables ~ ., scales = "free_y", space = "free_y")

ggsave(file.path(results_dir, "tissues_median_boxplot12.pdf"), plot = boxplot_mean, width = 12, height = plot_height)



# ------------------------------------------------------------------------------------



# Perform the t-test and save the p-values

t_test_results <- filtered_data %>% filter(!is.na(distance)) %>%
  group_by(tissue) %>%
  summarise(
    t_test = list(t.test(distance ~ type, data = ., alternative = 'greater')),
    .groups = "drop"  # To avoid warnings
  ) %>%
  ungroup() %>%
  mutate(
    p_value = map_dbl(t_test, ~ .x$p.value),        # Extract the p-value from the t-test
    statistic = map_dbl(t_test, ~ .x$statistic)     # Extract the test statistic from the t-test
  ) %>%
  select(-t_test)  # Remove the column with t-test objects to clean up results

# Apply p-value correction
t_test_results <- t_test_results %>%
  mutate(
    p_value_adjusted = p.adjust(p_value, method = "BH"), 
    significance_original = sapply(p_value, significance_level), 
    significance_adjusted = sapply(p_value_adjusted, significance_level)  
  )

write.csv(t_test_results, results_dir+"/t_test_results_by_tissue.csv", row.names = FALSE)

# Check that both groups are present
p <- ggplot(filtered_data, aes(x = tissue, y = distance, fill = type)) +
  geom_boxplot() +
  labs(title = "Expression Distances by Tissue",
       x = "Tissue",
       y = "Distance") +
  scale_fill_manual(values = c("Orthologs" = "blue", "Paralogs" = "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

for (i in 1:nrow(t_test_results)) {
  p <- p + annotate("text", 
                    x = i,  
                    y = max(filtered_data$distance, na.rm = TRUE) + 0.02, 
                    label = t_test_results$significance_adjusted[i], 
                    size = 5, 
                    color = "black")  
}

pdf(results_dir+"/boxplot_expression_distances_by_tissue.pdf", height = 15, width = 7)
print(p)
dev.off()





# ------------------------------------------------------------------------------------



# Perform Wilcoxon test for each tissue (variables) across gene types
stat.test <- mydata_long %>% 
  group_by(variables) %>% 
  wilcox_test(value ~ gene.type, p.adjust.method = "BH") %>% 
  add_significance()

# Create the boxplot
myplot <- ggboxplot(
  mydata_long,
  x = "gene.type",
  y = "value",
  fill = "gene.type",
  palette = "npg",  
  legend = "none",
  ggtheme = theme_pubr(border = TRUE)
) + facet_wrap(~variables)

# Add p-values to the plot
stat.test <- stat.test %>% add_xy_position(x = "gene.type")
myplot_with_pvals <- myplot + stat_pvalue_manual(stat.test, label = "p.adj.signif")

# Save the plot as a PDF
ggsave("plots/wilcox-test-plots/tissues_mean_boxplot_with_pvalues.pdf", plot = myplot_with_pvals, width = 12, height = 8)

# Save the statistical results as a .tsv file
write_tsv(stat.test, "plots/wilcox-test-plots/tissues_mean_statistical_results.tsv")