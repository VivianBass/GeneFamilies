
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
# Euclidean Distances regular (with mean/median statistics) 
# ------------------------------------------------------------------------

# Process regular distance statistics
regular_dfs <- process_distance_statistics(
    input_file = "exp.prof.dists_statistics.RData",
    stats_pattern = ".lst_dists_stats$",
    output_file = "exp.prof.dists_mean_median.RData",
    output_data_dir = output_data_dir
)

# Create plots for regular euclidean distance data
types <- unique(regular_dfs$mean$Type)
type_combinations <- combn(types, 2, simplify = FALSE)
plots_regular <- create_and_save_distance_boxplots(
    mean_data = regular_dfs$mean,
    median_data = regular_dfs$median,
    type_combinations = type_combinations,
    results_dir = results_dir,
    filename = "boxplots_euclidean_distances_regular.pdf"
)

# ------------------------------------------------------------------------
# Euclidean Distances log2 (with mean/median statistics)
# ------------------------------------------------------------------------

# Process log2 distance statistics
log2_dfs <- process_distance_statistics(
    input_file = "exp.prof.dists_statistics_log2.RData",
    stats_pattern = ".lst_dists_log2_stats$",
    output_file = "exp.prof.dists_mean_median_log2.RData",
    output_data_dir = output_data_dir
)

# Create plots for euclidean distance based on log2 transformed expression profils data
types <- unique(log2_dfs$mean$Type)
type_combinations <- combn(types, 2, simplify = FALSE)
plots_log2 <- create_and_save_distance_boxplots(
    mean_data = log2_dfs$mean,
    median_data = log2_dfs$median,
    type_combinations = type_combinations,
    results_dir = results_dir,
    filename = "boxplots_euclidean_distances_log2.pdf"
)

# ------------------------------------------------------------------------
# Euclidean Distances All (without mean/median statistics)
# ------------------------------------------------------------------------

# For regular data 
load(file.path(output_data_dir, "exp.prof.dists.RData"))
loaded_objects <- ls()
df_complete_dists <- process_regular_distances(".lst_dists$", loaded_objects)

# For log2 data
load(file.path(output_data_dir, "exp.prof.dists.log2.RData"))
loaded_objects <- ls()
df_complete_dists_log2 <- process_regular_distances(".lst_dists_log2$", loaded_objects)

combined_plots <- plot_and_save_combined_distances(
    df_complete_dists,
    df_complete_dists_log2,
    results_dir,
    "boxplots_euclidean_distances_all"
)







