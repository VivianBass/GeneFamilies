
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

message("USAGE: Rscript exec/plot_exp.prof.dists_euclidean_distribution.R")

# ------------------------------------------------------------------------
# Analysis and Visualization of Regular Euclidean Distances
# ------------------------------------------------------------------------

# Process statistics (mean and median) for regular Euclidean distances
regular_euclidean <- process_distances(
    input_file = "exp.prof.dists_statistics.RData",
    stats_pattern = ".lst_dists_stats$",
    output_file = "exp.prof.dists_mean_median.RData",
    output_data_dir = output_data_dir
)

# Perform statistical tests (t-test, Wilcoxon test) for computed distances
# Save the results in a CSV summary
regular_results <- perform_statistical_tests(
    regular_euclidean, 
    "statistical_tests_summary_euclidean.csv"
)

# Generate and save boxplots for Euclidean distance distributions
types <- unique(regular_euclidean$mean$Type) # Extract unique sample types
type_combinations <- combn(types, 2, simplify = FALSE) # Generate all pairwise comparisons
plots_regular <- create_distance_boxplots(
    mean_data = regular_euclidean$mean,
    median_data = regular_euclidean$median,
    type_combinations = type_combinations,
    metric = "Distance",
    results_dir = results_dir,
    filename = "boxplots_euclidean_distances_regular.pdf" # Output PDF file
)

# ------------------------------------------------------------------------
# Analysis and Visualization of Log2-Transformed Euclidean Distances
# ------------------------------------------------------------------------

# Process statistics (mean and median) for regular Euclidean distances
log2_euclidean <- process_distances(
    input_file = "exp.prof.dists_statistics_log2.RData",
    stats_pattern = ".lst_dists_log2_stats$",
    output_file = "exp.prof.dists_mean_median_log2.RData",
    output_data_dir = output_data_dir
)

# Perform statistical tests (t-test, Wilcoxon test) for log2-transformed distances
# Save the results in a CSV summary
log2_results <- perform_statistical_tests(
    log2_euclidean, 
    "statistical_tests_summary_euclidean_log2.csv"
)

# Generate and save boxplots for log2-transformed Euclidean distances
types <- unique(log2_euclidean$mean$Type) # Extract unique sample types
type_combinations <- combn(types, 2, simplify = FALSE) # Generate all pairwise comparisons
plots_log2 <- create_distance_boxplots(
    mean_data = log2_euclidean$mean,
    median_data = log2_euclidean$median,
    type_combinations = type_combinations,
    metric = "Distance",
    results_dir = results_dir,
    filename = "boxplots_euclidean_distances_log2.pdf" # Output PDF file
)

# ------------------------------------------------------------------------
# Analysis of Complete Euclidean Distances Without Aggregation (mean/median)
# ------------------------------------------------------------------------

# Process complete regular Euclidean distances without calculating mean/median
load(file.path(output_data_dir, "exp.prof.dists.RData")) # Load full data
loaded_objects <- ls() # List loaded objects
complete_dists <- process_complete_distances(".lst_dists$", loaded_objects)

# Perform statistical tests (t-test, Wilcoxon test) on complete data without aggregation
regular_results_all <- perform_statistical_tests_complete(
    complete_dists, 
    "statistical_tests_summary_euclidean_complete_(without_mean_median).csv"
)

# Process complete log2-transformed Euclidean distances without aggregation
load(file.path(output_data_dir, "exp.prof.dists.log2.RData")) # Load log2-transformed data
loaded_objects <- ls() # List loaded objects
complete_dists_log2 <- process_complete_distances(".lst_dists_log2$", loaded_objects)


# Perform statistical tests (t-test, Wilcoxon test) on log2-transformed data without aggregation
log2_results_all <- perform_statistical_tests_complete(
    complete_dists_log2, 
    "statistical_tests_summary_euclidean_complete_(without_mean_median)_log2.csv"
)

# Execute with combined data directly
combined_euclidean_plots <- create_distance_boxplots_all(
    regular_data = complete_dists,  
    log2_data = complete_dists_log2,  
    results_dir = results_dir,
    filename_prefix = "boxplots_euclidean_distances_all",
    metric = "Distance"
)



