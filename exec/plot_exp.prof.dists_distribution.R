require(GeneFamilies)
options(mc.cores = getMcCores())  
library(dotenv)  

# Load necessary libraries for data manipulation, plotting, and statistical analysis
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(rstatix)
library(ggpubr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

# Define directories for output data and results using environment variables
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

message("USAGE: Rscript exec/plot_exp.prof.dists_distribution.R")

# Load expression profile distance statistics
load(file.path(output_data_dir, "exp.prof.dists_statistics.RData"))
loaded_objects <- ls()
data_names <- loaded_objects[grepl(".lst_dists_stats$", loaded_objects)]

# Initialize empty data frames to store the mean and median distances for each dataset
df_mean.dists <- data.frame()
df_median.dists <- data.frame()

# Process each dataset to extract and compile mean and median distance information
for (data_name in data_names) {
    current_data <- get(data_name)
    type_name <- sub(".lst_dists_stats", "", data_name)
    
    # Compile mean distance data into a tibble, adding columns for type, cluster, and distance
    temp_mean_df <- tibble(
        Type = type_name,
        Cluster = names(current_data$Mean),
        Distance = unlist(current_data$Mean)
    )
    df_mean.dists <- bind_rows(df_mean.dists, temp_mean_df)
    
    # Compile median distance data into a tibble, adding columns for type, cluster, and distance
    temp_median_df <- tibble(
        Type = type_name,
        Cluster = names(current_data$Median),
        Distance = unlist(current_data$Median)
    )
    df_median.dists <- bind_rows(df_median.dists, temp_median_df)
}

# -----------------------------------------------------------------------------

# Function to label significance levels based on p-values
significance_level <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")  
}

# Check if there is enough data to perform t-tests on both mean and median distances
# Ensure each group has more than one observation for valid t-testing
valid_groups <- bind_rows(
    df_mean.dists %>% mutate(source = "mean"),
    df_median.dists %>% mutate(source = "median")
) %>%
    group_by(Type, source) %>%
    summarise(n = n(), .groups = 'drop') %>%
    filter(n > 1) %>%
    split(.$source) %>%
    map(~pull(.x, Type))

# Perform t-tests with error handling to account for potential issues
t_test_results <- tryCatch({

    # Perform t-test on median distances if there are enough groups
    if (length(valid_groups$median) >= 2) {
        t_test_median <- df_median.dists %>%
            filter(Type %in% valid_groups$median) %>%
            t_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(significance = sapply(p, significance_level),
                   analysis = "Median")
    } else {
        message("Not enough valid groups for median t-test")
    }
    
    # Perform t-test on mean distances if there are enough groups
    if (length(valid_groups$mean) >= 2) {
        t_test_mean <- df_mean.dists %>%
            filter(Type %in% valid_groups$mean) %>%
            t_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(significance = sapply(p, significance_level),
                   analysis = "Mean")
    } else {
        message("Not enough valid groups for mean t-test")
    }
    
    # Save results if both t-tests were completed successfully
    if (exists("t_test_median") && exists("t_test_mean")) {
        t_test_summary <- bind_rows(t_test_median, t_test_mean)
        write.csv(t_test_summary, file.path(results_dir, "t_test_summary.csv"), row.names = FALSE)
        message("T-test summary exported to CSV")
    }
    
    # Return results as a list for further reference
    list(median = if(exists("t_test_median")) t_test_median else NULL,
         mean = if(exists("t_test_mean")) t_test_mean else NULL,
         summary = if(exists("t_test_summary")) t_test_summary else NULL)
}, error = function(e) {
    message("Error in t-tests: ", e$message)
    return(NULL)
})

# -----------------------------------------------------------------------------

# Define combinations of types for pairwise significance annotations in boxplots
types <- unique(df_mean.dists$Type)
type_combinations <- combn(types, 2, simplify = FALSE)

# Plot the boxplot with significance annotations for mean distances
boxplot_mean <- ggplot(df_mean.dists, aes(x = Type, y = Distance, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Mean Expression Distances", y = "Distance") +
  theme_pubr(border = TRUE) +
  scale_y_continuous(breaks = seq(0, max(df_mean.dists$Distance), by = 0.2)) +
  geom_signif(comparisons = type_combinations, map_signif_level = TRUE) +
  theme(
    plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
    axis.text.x = element_text(size = 10),   
    plot.margin = margin(r = 30)
  )

# Plot the boxplot with significance annotations for median distances
boxplot_median <- ggplot(df_median.dists, aes(x = Type, y = Distance, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Median Expression Distances", y = "Distance") +
  theme_pubr(border = TRUE) +
  scale_y_continuous(breaks = seq(0, max(df_median.dists$Distance), by = 0.2)) +
  geom_signif(comparisons = type_combinations, map_signif_level = TRUE) +
  theme(
    plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
    axis.text.x = element_text(size = 10),   
    plot.margin = margin(r = 30)
  )

# Save both plots to a multi-page PDF
output_pdf <- file.path(results_dir, "boxplots_expression_distances.pdf")

# Combine plots into a list
plot_list <- list(boxplot_median, boxplot_mean)

# Save all plots with minimal page numbers
ggsave(output_pdf, marrangeGrob(plot_list, nrow=1, ncol=1, top=""), width = 12, height = 8)

