
message("USAGE: Rscript exec/generate_t-test_wilcox_test_tissue.R")

library(dotenv)
# Define directories for output data and results using environment variables
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Librarys for handling Dataframes, Lists etc. more efficiently
library(dplyr)
library(tidyr)
library(purrr)
library(parallel)

# Librarys for calculating t-tests and wilcox tests
library(rstatix)

# Functions sourced from:
source("R/compute_funks.R")

# ------------------------------------------------------------------------

# Load precomputed dataframes for mean and median distances generated in 
# the script `plot_exp.prof.dists_distribution.R`
load(file.path(output_data_dir, "exp.prof.dists_mean_median_tissue.RData"))

# Validate sufficient data for statistical tests
# Filter for groups with more than 1 observation per Tissue and Type 
# and organize by source (mean or median)
valid_groups_tissue <- bind_rows(
    df_mean.dists_tissue %>% mutate(source = "mean"),
    df_median.dists_tissue %>% mutate(source = "median")
) %>%
    group_by(Tissue, Type, source) %>%
    summarise(n = n(), .groups = 'drop') %>%
    filter(n > 1) %>%
    split(.$source) %>%
    map(~select(.x, Tissue, Type))

# ------------------------------------------------------------------------

# Perform statistical tests
# Conduct t-tests and Wilcoxon tests for mean and median distances to compare types 
# (e.g., orthologs vs. paralogs). Adjust p-values and summarize results.

test_results_tissue <- tryCatch({
    # Perform tests for median distances
    median_results <- perform_tissue_tests(df_median.dists_tissue, valid_groups_tissue, "median")
    # Perform tests for mean distances
    mean_results <- perform_tissue_tests(df_mean.dists_tissue, valid_groups_tissue, "mean")
    
    # Combine and save test summaries if results are available
    if (!is.null(median_results) && !is.null(mean_results)) {
        test_summary_tissue <- bind_rows(
            median_results$t_test, mean_results$t_test,
            median_results$wilcox, mean_results$wilcox
        )
        write.csv(test_summary_tissue, 
                  file.path(results_dir, "statistical_tests_summary_tissue1.csv"), 
                  row.names = FALSE)
        message("Tissue-specific statistical tests summary exported to CSV")
    }
    
    # Organize results for output
    list(
        t_test = list(
            median = if (!is.null(median_results)) median_results$t_test else NULL,
            mean = if (!is.null(mean_results)) mean_results$t_test else NULL
        ),
        wilcox = list(
            median = if (!is.null(median_results)) median_results$wilcox else NULL,
            mean = if (!is.null(mean_results)) mean_results$wilcox else NULL
        ),
        summary = if (exists("test_summary_tissue")) test_summary_tissue else NULL
    )
}, error = function(e) {
    message("Error in tissue-specific statistical tests: ", e$message)
    return(NULL)
})
