
# Load environment variables to define directories for output data and results
library(dotenv)
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Create results directory if it doesn't exist
if(!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}

message("USAGE: Rscript exec/plot_exp.prof.dists_distribution.R")

# Libraries for efficient data handling
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(parallel)

# Libraries for data visualization
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(ggpubr)

# functions sourced from:
source("R/plot_distribution_funks.R")

# ------------------------------------------------------------------------

# Process regular data
regular_data <- process_tissue_statistics(
    file.path(output_data_dir, "exp.prof.dists_statistics.RData"),
    "_tissue_stats$",
    ".lst_dists_tissue_stats"
)

# Process log2 data
log2_data <- process_tissue_statistics(
    file.path(output_data_dir, "exp.prof.dists_statistics_log2.RData"),
    "_tissue_log2_stats$",
    ".lst_dists_tissue_log2_stats"
)

# Save processed data
save(regular_data$mean, regular_data$median,
     file = file.path(output_data_dir, "exp.prof.dists_mean_median_tissue.RData"))
save(log2_data$mean, log2_data$median,
     file = file.path(output_data_dir, "exp.prof.dists_mean_median_tissue_log2.RData"))

# Generate all plots
for (test_type in c("t.test", "wilcox.test")) {
    for (metric_type in c("Mean", "Median")) {
        # Regular data plots
        create_tissue_boxplots(
            if(metric_type == "Mean") regular_data$mean else regular_data$median,
            metric_type, test_type, FALSE, results_dir
        )
        
        # Log2 data plots
        create_tissue_boxplots(
            if(metric_type == "Mean") log2_data$mean else log2_data$median,
            metric_type, test_type, TRUE, results_dir
        )
    }
}

# ------------------------------------------------------------------------
# Euclidean Distances All (without mean/median statistics)
# ------------------------------------------------------------------------

# For regular data 
load(file.path(output_data_dir, "exp.prof.dists.RData"))
loaded_objects <- ls()
df_complete_dists_tissue <- process_tissue_distances(".lst_dists_tissue$", loaded_objects)

# For log2 data
load(file.path(output_data_dir, "exp.prof.dists.log2.RData"))
loaded_objects <- ls()
df_complete_dists_tissue_log2 <- process_tissue_distances(".lst_dists_tissue_log2$", loaded_objects)


plot_and_save_combined_tissue_distances(
    df_complete_dists_tissue,
    df_complete_dists_tissue_log2,
    results_dir,
    "tissues_all_boxplot"
)

# Generate plots for regular data
create_tissue_boxplots(df_complete_dists_tissue, "All", "t.test", FALSE, results_dir)
create_tissue_boxplots(df_complete_dists_tissue, "All", "wilcox.test", FALSE, results_dir)

# Generate plots for log2 data
create_tissue_boxplots(df_complete_dists_tissue_log2, "All", "t.test", TRUE, results_dir)
create_tissue_boxplots(df_complete_dists_tissue_log2, "All", "wilcox.test", TRUE, results_dir)
