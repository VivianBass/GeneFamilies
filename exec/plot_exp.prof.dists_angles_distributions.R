
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
# Angles Distances regular (with mean/median statistics) 
# ------------------------------------------------------------------------

# Process angles distance statistics
angles_dfs <- process_distance_statistics(
    input_file = "exp.prof.angels_statistics.RData",
    stats_pattern = ".lst_cos_angles_dists_stats$",
    output_file = "exp.prof.angles_mean_median.RData",
    output_data_dir = output_data_dir
)

# Create plots for cosine angles distances data
types <- unique(angles_dfs$mean$Type)
type_combinations <- combn(types, 2, simplify = FALSE)
plots_angles <- create_and_save_distance_boxplots(
    mean_data = angles_dfs$mean,
    median_data = angles_dfs$median,
    type_combinations = type_combinations,
    results_dir = results_dir,
    filename = "boxplots_angles_distances_regular.pdf"
)

# ------------------------------------------------------------------------
# Angles Distances log2 (with mean/median statistics) 
# ------------------------------------------------------------------------

# Process angles log2 distance statistics
angles_log2_dfs <- process_distance_statistics(
    input_file = "exp.prof.angels_statistics_log2.RData",
    stats_pattern = ".lst_cos_angles_dists_log2_stats$",
    output_file = "exp.prof.angles_mean_median_log2.RData",
    output_data_dir = output_data_dir
)

# Create plots for log2 transformed cosine angles distances data
types <- unique(angles_log2_dfs$mean$Type)
type_combinations <- combn(types, 2, simplify = FALSE)
plots_angles_log2 <- create_and_save_distance_boxplots(
    mean_data = angles_log2_dfs$mean,
    median_data = angles_log2_dfs$median,
    type_combinations = type_combinations,
    results_dir = results_dir,
    filename = "boxplots_angles_distances_log2.pdf"
)

# ------------------------------------------------------------------------
# Angles Distances All (without mean/median statistics)
# ------------------------------------------------------------------------

# Usage:
# For regular angles
load(file.path(output_data_dir, "exp.prof.angles.RData"))
loaded_objects <- ls()
df_complete_angles <- process_regular_angles(".lst_cos_angles_dists$", loaded_objects)

# For log2 angles
load(file.path(output_data_dir, "exp.prof.angles.log2.RData"))
loaded_objects <- ls()
df_complete_angles_log2 <- process_regular_angles(".lst_cos_angles_dists_log2$", loaded_objects)

combined_angle_plots <- plot_and_save_combined_angles(
    df_complete_angles,
    df_complete_angles_log2,
    results_dir,
    "boxplots_angles_distances_all"
)