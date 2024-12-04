
# Import necessary libraries
library(dotenv) # Load environment variables for directory paths
library(dplyr) # Data manipulation
library(tidyr) # Data tidying
library(purrr) # Functional programming for data
library(tibble) # Handling tidy data frames
library(parallel) # Parallel processing for performance optimization
library(ggplot2) # Data visualization
library(ggsignif) # Add significance annotations to ggplots
library(gridExtra) # Arrange multiple plots
library(ggpubr) # Publication-ready ggplots
library(rstatix) # Perform statistical tests (t-tests, Wilcoxon tests)

# Load user-defined functions
source("R/compute_funks.R") # Functions for computing distance statistics
source("R/statistical_tests_funks.R") # Functions for performing statistical tests
source("R/plot_distribution_funks.R") # Functions for creating and saving plots

# Set output directories from environment variables
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Ensure the results directory exists, create it if not
if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}

message("USAGE: Rscript exec/plot_exp.prof.dists_angles_distribution.R")

# ------------------------------------------------------------------------
# Angles Distances regular (with mean/median statistics) 
# ------------------------------------------------------------------------

# Process angles distance statistics
regular_angles <- process_distances(
    input_file = "exp.prof.angels_statistics.RData",
    stats_pattern = ".lst_cos_angles_dists_stats$",
    output_file = "exp.prof.angles_mean_median.RData",
    output_data_dir = output_data_dir
)

# Perform statistical tests (t-test, Wilcoxon test) for computed distances
# Save the results in a CSV summary
regular_results <- perform_statistical_tests(regular_angles, "statistical_tests_summary_angles.csv")

# Generate type combinations before creating plots
types <- unique(regular_angles$mean$Type)
type_combinations <- combn(types, 2, simplify = FALSE)
plots_angles <- create_distance_boxplots(
    mean_data = regular_angles$mean,
    median_data = regular_angles$median,
    type_combinations = type_combinations,
    metric = "Angle",
    results_dir = results_dir,
    filename = "boxplots_angles_distances_regular.pdf"
)

# ------------------------------------------------------------------------
# Angles Distances log2 (with mean/median statistics) 
# ------------------------------------------------------------------------

# Process angles log2 distance statistics
log2_angles <- process_distances(
    input_file = "exp.prof.angels_statistics_log2.RData",
    stats_pattern = ".lst_cos_angles_dists_log2_stats$",
    output_file = "exp.prof.angles_mean_median_log2.RData",
    output_data_dir = output_data_dir
)

regular_results <- perform_statistical_tests(log2_angles, "statistical_tests_summary_angles_log2.csv")

# Create plots for log2 transformed cosine angles distances data
types <- unique(log2_angles$mean$Type)
type_combinations <- combn(types, 2, simplify = FALSE)
plots_angles_log2 <- create_distance_boxplots(
    mean_data = log2_angles$mean,
    median_data = log2_angles$median,
    type_combinations = type_combinations,
    metric = "Angle",
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
complete_angles <- process_complete_distances(".lst_cos_angles_dists$", loaded_objects)

# Perform statistical tests (t-test, Wilcoxon test) on complete data without aggregation
complete_angles_results <- perform_statistical_tests_complete(
    complete_angles, 
    "statistical_tests_summary_angles_complete_(without_mean_median).csv"
)

# For log2 angles
load(file.path(output_data_dir, "exp.prof.angles.log2.RData"))
loaded_objects <- ls()
complete_log2_angles <- process_complete_distances(".lst_cos_angles_dists_log2$", loaded_objects)

# Perform statistical tests (t-test, Wilcoxon test) on complete data without aggregation
complete_log2_angles_results  <- perform_statistical_tests_complete(
    complete_log2_angles, 
    "statistical_tests_summary_angles_log2_complete_(without_mean_median).csv"
)

# Execute with combined data directly
combined_angle_plots <- create_distance_boxplots_all(
    regular_data = complete_angles,  # Complete dataset
    log2_data = complete_log2_angles,  # Complete log2 dataset
    results_dir = results_dir,
    filename_prefix = "boxplots_angles_distances_all",
    metric = "Angle"
)









# ------------------------------------------------------------------------
# Angles Distances regular in degrees (with mean/median statistics) 
# ------------------------------------------------------------------------

# Process angles distance statistics for degree data
regular_angles_degree <- process_distances(
    input_file = "exp.prof.angels_statistics_degree.RData",
    stats_pattern = ".lst_cos_angles_dists_degree_stats$",
    output_file = "exp.prof.angles_mean_median_degree.RData",
    output_data_dir = output_data_dir
)

# Perform statistical tests (t-test, Wilcoxon test) for computed distances in degrees
# Save the results in a CSV summary
regular_results_degree <- perform_statistical_tests(regular_angles_degree, "statistical_tests_summary_angles_degree.csv")

# Generate type combinations before creating plots for degree data
types_degree <- unique(regular_angles_degree$mean$Type)
type_combinations_degree <- combn(types_degree, 2, simplify = FALSE)
plots_angles_degree <- create_distance_boxplots(
    mean_data = regular_angles_degree$mean,
    median_data = regular_angles_degree$median,
    type_combinations = type_combinations_degree,
    metric = "Angle (Degrees)",
    results_dir = results_dir,
    filename = "boxplots_angles_distances_regular_degree.pdf"
)

# ------------------------------------------------------------------------
# Angles Distances log2 in degrees (with mean/median statistics) 
# ------------------------------------------------------------------------

# Process angles log2 distance statistics for degree data
log2_angles_degree <- process_distances(
    input_file = "exp.prof.angels_statistics_log2_degree.RData",
    stats_pattern = ".lst_cos_angles_dists_log2_degree_stats$",
    output_file = "exp.prof.angles_mean_median_log2_degree.RData",
    output_data_dir = output_data_dir
)

# Perform statistical tests for log2 transformed degree data
log2_results_degree <- perform_statistical_tests(log2_angles_degree, "statistical_tests_summary_angles_log2_degree.csv")

# Create plots for log2 transformed cosine angles distances data in degrees
types_log2_degree <- unique(log2_angles_degree$mean$Type)
type_combinations_log2_degree <- combn(types_log2_degree, 2, simplify = FALSE)
plots_angles_log2_degree <- create_distance_boxplots(
    mean_data = log2_angles_degree$mean,
    median_data = log2_angles_degree$median,
    type_combinations = type_combinations_log2_degree,
    metric = "Angle (Degrees)",
    results_dir = results_dir,
    filename = "boxplots_angles_distances_log2_degree.pdf"
)

