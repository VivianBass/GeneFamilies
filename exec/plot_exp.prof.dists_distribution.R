require(GeneFamilies)
options(mc.cores = getMcCores())

library(ggplot2)
library(ggsignif)
library(gridExtra)
library(dotenv)

library(tidyr)
library(dplyr)
library(purrr)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")


load(file.path(output_data_dir, "expression_profile_distances_statistics.RData"))
load(file.path(output_data_dir, "ExpressionProfileDistances.RData"))

# Function to determine the level of significance
significance_level <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")  
}

# - combined df ???? from genefamilies

# mean and median expression distances
# df_median_mean_paralogs 
# df_median_mean_orthologs

paralog.mean.lst <- split(df_median_mean_paralogs$Mean, df_median_mean_paralogs$Family)
paralog.median.lst <- split(df_median_mean_paralogs$Median, df_median_mean_paralogs$Family)

ortholog.mean.lst <- split(df_median_mean_orthologs$Mean, df_median_mean_orthologs$Family)
ortholog.median.lst <- split(df_median_mean_orthologs$Median, df_median_mean_orthologs$Family)


# basically create df from the lists, beneath each other in the df, regardless of df length
df_mean.dists <- map_dfr(list(Ortholog = ortholog.mean.lst, Paralog = paralog.mean.lst), 
              ~tibble(Cluster = names(.x), Distance = unlist(.x)), .id = "Type")

df_median.dists <- map_dfr(list(Ortholog = ortholog.median.lst, Paralog = paralog.median.lst), 
              ~tibble(Cluster = names(.x), Distance = unlist(.x)), .id = "Type")


# Create a boxplot for mean distances
boxplot_mean <- ggplot(df_mean.dists, aes(x = Type, y = Distance, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Mean Expression Distances", y = "Distance") +
  theme_minimal() +
  geom_signif(comparisons = list(c("Ortholog", "Paralog")), map_signif_level = TRUE) +
  theme(plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, r = 0 ,b = 20, l = 0), hjust = 0.50),
        axis.title.x = element_text(size = 10, margin = margin(t = 20, r = 0, b = 20, l = 0), hjust = 0.50),
        axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
  )

# Create a boxplot for median distances
boxplot_median <- ggplot(df_median.dists, aes(x = Type, y = Distance, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +  # Exclude outliers from the boxplot to avoid overlap
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +  # Add jittered points
  labs(title = "Median Expression Distances", y = "Distance") +
  theme_minimal() +
  geom_signif(comparisons = list(c("Orthologs", "Paralogs")), map_signif_level = TRUE) +
  theme(plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, r = 0 ,b = 20, l = 0), hjust = 0.50),
        axis.title.x = element_text(size = 10, margin = margin(t = 20, r = 0, b = 20, l = 0), hjust = 0.50),
        axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
  )

# Save all three boxplots in a single PDF
pdf(file.path(results_dir,"boxplots_expression_distances_with_jitter.pdf"), height = 15, width = 7)
grid.arrange(boxplot_median, boxplot_mean, ncol = 1)
dev.off()

# ----------------------------------------------------------------

# assuming they have the same row length ?? i guess ??
# assuming gene_family_expression_dists_stats is a combined df ??
# should we use this code ??

# Create a dataframe with median distances
median_distances <- data.frame(
  distance = c(gene_family_expression_dists_stats$median_ortholog_exp_dists,
               gene_family_expression_dists_stats$median_paralog_exp_dists),
  group = rep(c("Orthologs", "Paralogs"), each = nrow(gene_family_expression_dists_stats))
)

# Create a dataframe with mean distances
mean_distances <- data.frame(
  distance = c(gene_family_expression_dists_stats$mean_ortholog_exp_dists,
               gene_family_expression_dists_stats$mean_paralog_exp_dists),
  group = rep(c("Orthologs", "Paralogs"), each = nrow(gene_family_expression_dists_stats))
)

# t-tests
# # Uncomment this part to perform t-test for medians
# t_test_medians <- t.test(distance ~ group, data = median_distances, alternative = 'greater')

# Uncomment this part to perform t-test for means
# t_test_mean_results <- t.test(distance ~ group, data = mean_distances, alternative = 'greater')

# Uncomment this part to save the t-test results in a dataframe
# t_test_summary <- data.frame(
#   analysis = c("Median", "Mean"),

  
#   # Results without correction
#   statistic_original = c(
#                           t_test_medians$statistic, 
#                           t_test_mean_results$statistic),
#   p_value_original = c(
#                        t_test_medians$p.value, 
#                        t_test_mean_results$p.value),
#   conf_int_lower_original = c(
#                                t_test_medians$conf.int[1], 
#                                t_test_mean_results$conf.int[1]),
#   conf_int_upper_original = c(
#                                t_test_medians$conf.int[2], 
#                                t_test_mean_results$conf.int[2]),
#   # Use the means reported by the t-tests
#   mean_orthologs = c(
#                       t_test_medians$estimate[1],     
#                       t_test_mean_results$estimate[1]),  
#   mean_paralogs = c(
#                      t_test_medians$estimate[2],     
#                      t_test_mean_results$estimate[2]), 
  
#   # Significance
#   significance_original = c(
#                     significance_level(t_test_medians$p.value),
#                     significance_level(t_test_mean_results$p.value)),
#   p_value_corrected = p.adjust(c(
#                                   t_test_medians$p.value, 
#                                   t_test_mean_results$p.value), 
#                                 method = "BH")
# )

# # Add a column of significance based on the corrected p-value
# t_test_summary$significance_corrected <- sapply(t_test_summary$p_value_corrected, significance_level)

# # Uncomment this part to save the original and corrected t-test results to a CSV file
# write.csv(t_test_summary, file.path(results_dir, "t_test_results_original_and_corrected.csv"), row.names = FALSE)




