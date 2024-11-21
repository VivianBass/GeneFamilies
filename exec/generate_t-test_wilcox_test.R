require(GeneFamilies)
options(mc.cores = getMcCores())
library(parallel)

message("USAGE: Rscript exec/generate_t-test_wilcox_test.R")

library(dotenv)
# Define directories for output data and results using environment variables
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Librarys for handling Dataframes, Lists etc. more efficiently
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

# Librarys for calculating t-tests and wilcox tests
library(rstatix)

# ------------------------------------------------------------------------

# load the 2 dataframes for mean and median distances geneerated in 
# exec/plot_exp.prof.dists_distribution.R 
load(file.path(output_data_dir, "exp.prof.dists_mean_median.RData"))

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

# ------------------------------------------------------------------------

# Function to label significance levels based on p-values
significance_level <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")  
}

# ------------------------------------------------------------------------

# Perform both t-tests and Wilcoxon tests with error handling
test_results <- tryCatch({
    
    # Median tests
    if (length(valid_groups$median) >= 2) {
        t_test_median <- df_median.dists %>%
            filter(Type %in% valid_groups$median) %>%
            t_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(significance = sapply(p, significance_level),
                   analysis = "Median",
                   test_type = "t-test")
                   
        wilcox_median <- df_median.dists %>%
            filter(Type %in% valid_groups$median) %>%
            wilcox_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(significance = sapply(p, significance_level),
                   analysis = "Median",
                   test_type = "wilcox")
    }
    
    # Mean tests
    if (length(valid_groups$mean) >= 2) {
        t_test_mean <- df_mean.dists %>%
            filter(Type %in% valid_groups$mean) %>%
            t_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(significance = sapply(p, significance_level),
                   analysis = "Mean",
                   test_type = "t-test")
                   
        wilcox_mean <- df_mean.dists %>%
            filter(Type %in% valid_groups$mean) %>%
            wilcox_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(significance = sapply(p, significance_level),
                   analysis = "Mean",
                   test_type = "wilcox")
    }
    
    # Combine and save all test results
    if (exists("t_test_median") && exists("t_test_mean") && 
        exists("wilcox_median") && exists("wilcox_mean")) {
        test_summary <- bind_rows(
            t_test_median, t_test_mean,
            wilcox_median, wilcox_mean
        )
        write.csv(test_summary, file.path(results_dir, "statistical_tests_summary.csv"), row.names = FALSE)
        message("Statistical tests summary exported to CSV")
    }
    
    # Return results
    list(
        t_test = list(
            median = if(exists("t_test_median")) t_test_median else NULL,
            mean = if(exists("t_test_mean")) t_test_mean else NULL
        ),
        wilcox = list(
            median = if(exists("wilcox_median")) wilcox_median else NULL,
            mean = if(exists("wilcox_mean")) wilcox_mean else NULL
        ),
        summary = if(exists("test_summary")) test_summary else NULL
    )
}, error = function(e) {
    message("Error in statistical tests: ", e$message)
    return(NULL)
})
