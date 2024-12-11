
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
source("R/statistical_tests_tissue_funks.R") # Functions for performing statistical tests
source("R/plot_distribution_tissue_funks.R") # Functions for creating and saving plots

# Set output directories from environment variables
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Ensure the results directory exists, create it if not
if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}

message("USAGE: Rscript exec/plot_exp.prof.dists_euclidean_tissue_distributions.R")

# ------------------------------------------------------------------------
# Process and Analyze Euclidean Distance Statistics for Tissues
# ------------------------------------------------------------------------

# Process statistics (mean, median) for regular data
regular_data <- process_tissue_statistics(
    file.path(output_data_dir, "exp.prof.dists_statistics.RData"),
    "_tissue_stats$",
    ".lst_dists_tissue_stats"
)

# Process summary statistics (mean, median) for log2-transformed data
log2_data <- process_tissue_statistics(
    file.path(output_data_dir, "exp.prof.dists_statistics_log2.RData"),
    "_tissue_log2_stats$",
    ".lst_dists_tissue_log2_stats"
)

# Perform statistical tests (t-test, Wilcoxon test) on statistics (mean and median)
# and save results in a CSV summary

# For regular data
mean_data <- regular_data$mean
median_data <- regular_data$median

test_results_regular <- perform_tissue_statistical_analysis(mean_data, median_data, results_dir, is_log2 = FALSE)

# For log2 transformed data
log2_mean_data <- log2_data$mean
log2_median_data <- log2_data$median

test_results_log2 <- perform_tissue_statistical_analysis(log2_mean_data, log2_median_data, results_dir, is_log2 = TRUE)


# Generate and save boxplots comparing Euclidean distance distributions for tissues
create_tissue_boxplots_combined(
    mean_data = regular_data$mean,
    median_data = regular_data$median,
    log2_mean_data = log2_data$mean,
    log2_median_data = log2_data$median,
    results_dir = results_dir
)

# ------------------------------------------------------------------------
# Euclidean Distances All (without mean/median statistics)
# ------------------------------------------------------------------------

# Load and process full Euclidean distance for regular data
load(file.path(output_data_dir, "exp.prof.dists.RData"))
loaded_objects <- ls()
complete_dists_tissue <- process_tissue_distances(".lst_dists_tissue$", loaded_objects)

# Load and process full Euclidean distance for log2-transformed data
load(file.path(output_data_dir, "exp.prof.dists.log2.RData"))
loaded_objects <- ls()
complete_dists_tissue_log2 <- process_tissue_distances(".lst_dists_tissue_log2$", loaded_objects)


# Perform statistical tests (t-test, Wilcoxon test) on full Euclidean distance distributions
# and save results in a CSV summary

# For regular data
test_results_regular_complete <- perform_tissue_statistical_analysis_complete(
    complete_dists_tissue,
    "statistical_tests_summary_tissue_complete_regular.csv",
    results_dir,
    is_log2 = FALSE
)

# For log2 data
test_results_log2_complete <- perform_tissue_statistical_analysis_complete(
    complete_dists_tissue_log2,
    "statistical_tests_summary_tissue_complete_log2.csv", 
    results_dir,
    is_log2 = TRUE
)

# Generate and save boxplots for full Euclidean distance distributions
create_tissue_boxplots_all_combined(
    complete_dists_tissue,
    complete_dists_tissue_log2,
    results_dir, "tissue_boxplot_all"
)









