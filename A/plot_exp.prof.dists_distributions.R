
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
library(parallel)

# Libraries for data visualization
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(ggpubr)

# functions sourced from:
source("R/plot_distribution_funks.R")


# plot_angles_to_diagonal_distributions
- exec\plot_exp.prof.dists_angles.R


# plot_exp.prof.dists_distributions_all_(without_statistics) 
- exec\plot_exp.prof.angles_distribution_(all).R
- exec\plot_exp.prof.dists_distribution_(all).R
- exec\plot_exp.prof.dists_distribution_log2_(all).R
- exec\plot_exp.prof.angles_distribution_log2_(all).R


- functions -> R/plot_tissue_distributions_funks.R
# plot_exp.prof.dists_tissue_distributions 
- exec\plot_exp.prof.dists_distribution_tissue_log2.R
- exec\plot_exp.prof.dists_distribution_tissue.R






# ------------------------------------------------------------------------
# regular and log 2
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
    filename = "boxplots_expression_distances_regular_all_tests.pdf"
)

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
    filename = "boxplots_expression_distances_log2_all_tests.pdf"
)

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
    filename = "boxplots_expression_distances_angles_all_tests.pdf"
)

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
    filename = "boxplots_expression_distances_angles_log2_all_tests.pdf"
)

# ---------------------------------------------------------------------
# Tissue
# ---------------------------------------------------------------------

# For regular data 
load(file.path(output_data_dir, "exp.prof.dists.RData"))
loaded_objects <- ls()

df_complete_dists <- process_regular_distances(".lst_dists$", loaded_objects)
df_complete_dists_tissue <- process_tissue_distances(".lst_dists_tissue$", loaded_objects)

# For log2 data
load(file.path(output_data_dir, "exp.prof.dists.log2.RData"))
loaded_objects <- ls()

df_complete_dists_log2 <- process_regular_distances(".lst_dists_log2$", loaded_objects)
df_complete_dists_tissue_log2 <- process_tissue_distances(".lst_dists_tissue_log2$", loaded_objects)


combined_plots <- plot_and_save_combined_distances(
    df_complete_dists,
    df_complete_dists_log2,
    results_dir,
    "complete_boxplots_expression_distances"
)

plot_and_save_combined_tissue_distances(
    df_complete_dists_tissue,
    df_complete_dists_tissue_log2,
    results_dir,
    "complete_boxplots_tissues_expression_distances"
)
















# ---------------------------------------------------------------------
# All 
# ---------------------------------------------------------------------


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
    "complete_boxplots_expression_angles"
)



# ---------------------------------------------------------------------------

# For regular data 
load(file.path(output_data_dir, "exp.prof.dists.RData"))
loaded_objects <- ls()

df_complete_dists <- process_regular_distances(".lst_dists$", loaded_objects)
df_complete_dists_tissue <- process_tissue_distances(".lst_dists_tissue$", loaded_objects)

# For log2 data
load(file.path(output_data_dir, "exp.prof.dists.log2.RData"))
loaded_objects <- ls()

df_complete_dists_log2 <- process_regular_distances(".lst_dists_log2$", loaded_objects)
df_complete_dists_tissue_log2 <- process_tissue_distances(".lst_dists_tissue_log2$", loaded_objects)


combined_plots <- plot_and_save_combined_distances(
    df_complete_dists,
    df_complete_dists_log2,
    results_dir,
    "complete_boxplots_expression_distances"
)

plot_and_save_combined_tissue_distances(
    df_complete_dists_tissue,
    df_complete_dists_tissue_log2,
    results_dir,
    "complete_boxplots_tissues_expression_distances"
)

# ---------------------------------------------------------------------------