require(GeneFamilies)
options(mc.cores = getMcCores())
library(dotenv)

library(ggplot2)
library(ggsignif)
library(gridExtra)
library(dotenv)

library(tidyr)
library(dplyr)
library(purrr)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

load(file.path(output_data_dir, "exp.prof.dists_statistics.RData"))

# Track created dataframe names for plotting
# 5 gene groups , 1 ortholog and 4 paralogs, families, 
# [23] "special_in_paralogs_filtered.filtered.dists.tissue_stats"
# [24] "special_in_paralogs_filtered.filtered.dists_stats"

# Automatically sort the loadedstatistics data into regular and tissue datasets
load(file.path(output_data_dir, "exp.prof.dists_statistics.RData"))
loaded_objects <- ls()
data_names <- loaded_objects[grepl("\\.dists_stats$", loaded_objects)]
data_names_tissue <- loaded_objects[grepl("\\.dists\\.tissue_stats$", loaded_objects)]

# Initialize empty dataframes
df_mean.dists <- data.frame()
df_median.dists <- data.frame()

# Process each dataset
for (data_name in data_names) {

    current_data <- get(data_name)
    type_name <- sub("_filtered.filtered.dists_stats", "", data_name)
    
    # Create and combine mean data
    temp_mean_df <- tibble(
        Type = type_name,
        Cluster = names(current_data$Mean),
        Distance = unlist(current_data$Mean)
    )
    df_mean.dists <- bind_rows(df_mean.dists, temp_mean_df)
    
    # Create and combine median data
    temp_median_df <- tibble(
        Type = type_name,
        Cluster = names(current_data$Median),
        Distance = unlist(current_data$Median)
    )
    df_median.dists <- bind_rows(df_median.dists, temp_median_df)
}

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
# Arrange plots using gridExtra
combined_plots <- grid.arrange(boxplot_median, boxplot_mean, ncol = 1)
ggsave(filename = file.path(results_dir, "boxplots_expression_distances.pdf"),
       plot = combined_plots, height = 15, width = 7)





