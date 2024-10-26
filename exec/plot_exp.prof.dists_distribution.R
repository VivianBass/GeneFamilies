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

# mean and median expression distances
# df_median_mean_paralogs 
# df_median_mean_orthologs

paralog.mean.lst <- split(df_median_mean_paralogs$Mean, df_median_mean_paralogs$Family)
paralog.median.lst <- split(df_median_mean_paralogs$Median, df_median_mean_paralogs$Family)

ortholog.mean.lst <- split(df_median_mean_orthologs$Mean, df_median_mean_orthologs$Family)
ortholog.median.lst <- split(df_median_mean_orthologs$Median, df_median_mean_orthologs$Family)


# basically create df from the lists, beneath each other in the df, regardless of df length
df_mean.dists <- map_dfr(list(Ortholog = ortholog.mean.lst, Paralog = paralog.mean.lst), ~tibble(Cluster = names(.x), Distance = unlist(.x)), .id = "Type")

df_median.dists <- map_dfr(list(Ortholog = ortholog.median.lst, Paralog = paralog.median.lst), ~tibble(Cluster = names(.x), Distance = unlist(.x)), .id = "Type")


# ---------------------------------------------------------------------------

# t-tests pairwise distribution data in long format (already in long format)

# Function to determine the level of significance
significance_level <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")  
}

# Perform t-test for median and adjust p-values, then apply the significance_level function

t_test_median <- df_median.dists %>%
  t_test(Distance ~ Type, alternative = "greater") %>% 
  adjust_pvalue(method = "BH") %>%                    
  mutate(significance = sapply(p, significance_level),
         analysis = "Median")  

t_test_mean <- df_median.dists %>%
  t_test(Distance ~ Type, alternative = "greater") %>% 
  adjust_pvalue(method = "BH") %>%                    
  mutate(significance = sapply(p, significance_level),
         analysis = "Mean")  

# Combine the two results into one summary dataframe
t_test_summary <- bind_rows(t_test_median, t_test_mean)

write.csv(t_test_summary, file.path(results_dir, "t_test_summary.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------

# Extract the significance levels are already in the boxplot


# Create the boxplot with the significance annotation
boxplot_mean <- ggplot(df_mean.dists, aes(x = Type, y = Distance, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Mean Expression Distances", y = "Distance") +
  theme_minimal() +
  geom_signif(comparisons = list(c("Ortholog", "Paralog")), map_signif_level = TRUE) +
  theme(
    plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, r = 0 ,b = 20, l = 0), hjust = 0.50),
    axis.title.x = element_text(size = 10, margin = margin(t = 20, r = 0, b = 20, l = 0), hjust = 0.50),
    axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20))
  )

# Create the boxplot with the significance annotation for median distances
boxplot_median <- ggplot(df_median.dists, aes(x = Type, y = Distance, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Median Expression Distances", y = "Distance") +
  theme_minimal() +
  geom_signif(comparisons = list(c("Ortholog", "Paralog")), map_signif_level = TRUE) +
  theme(
    plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, r = 0 ,b = 20, l = 0), hjust = 0.50),
    axis.title.x = element_text(size = 10, margin = margin(t = 20, r = 0, b = 20, l = 0), hjust = 0.50),
    axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20))
  )

# Save the two Boxplots in a single PDF
pdf(file.path(results_dir,"boxplots_expression_distances_with_jitter11.pdf"), height = 15, width = 7)
grid.arrange(boxplot_median, boxplot_mean, ncol = 1)
dev.off()





