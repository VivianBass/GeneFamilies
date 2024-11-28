
library(dotenv)
# Define directories for output data and results using environment variables
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Create results directory if it doesn't exist
if(!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}

message("USAGE: Rscript exec/generate_t-test_wilcox_test_cos_angles_log2.R")

# Librarys for handling Dataframes, Lists etc. more efficiently
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(parallel)

# Librarys for calculating t-tests and wilcox tests
library(rstatix)

# Functions sourced from:
source("R/compute_funks.R")

# ------------------------------------------------------------------------

# load the 2 dataframes for mean and median distances
load(file.path(output_data_dir, "exp.prof.angles_mean_median_log2.RData"))

# Check if there is enough data to perform t-tests
valid_groups <- bind_rows(
    df_mean.dists %>% mutate(source = "mean"),
    df_median.dists %>% mutate(source = "median")
) %>%
    group_by(Type, source) %>%
    summarise(n = n(), .groups = 'drop') %>%
    filter(n > 1) %>%
    split(.$source) %>%
    map(~pull(.x, Type))

# ------------------------------------------------------------------------

# Perform tests with error handling
test_results <- tryCatch({
    median_results <- perform_tests(df_median.dists, valid_groups, "median")
    mean_results <- perform_tests(df_mean.dists, valid_groups, "mean")
    
    if (!is.null(median_results) && !is.null(mean_results)) {
        test_summary <- bind_rows(
            median_results$t_test, mean_results$t_test,
            median_results$wilcox, mean_results$wilcox
        )
        write.csv(test_summary, 
                 file.path(results_dir, "statistical_tests_summary_angles_log2.csv"), 
                 row.names = FALSE)
        message("Statistical tests summary for angles (log2) exported to CSV")
    }
    
    list(
        t_test = list(
            median = if(!is.null(median_results)) median_results$t_test else NULL,
            mean = if(!is.null(mean_results)) mean_results$t_test else NULL
        ),
        wilcox = list(
            median = if(!is.null(median_results)) median_results$wilcox else NULL,
            mean = if(!is.null(mean_results)) mean_results$wilcox else NULL
        ),
        summary = if(exists("test_summary")) test_summary else NULL
    )
}, error = function(e) {
    message("Error in statistical tests: ", e$message)
    return(NULL)
})



