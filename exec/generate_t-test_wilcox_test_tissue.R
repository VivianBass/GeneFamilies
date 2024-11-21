require(GeneFamilies)
options(mc.cores = getMcCores())
library(parallel)

message("USAGE: Rscript exec/generate_t-test_wilcox_test_tissue.R")

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
load(file.path(output_data_dir, "exp.prof.dists_mean_median_tissue.RData"))

# Check if sufficient data exists for t-tests by counting observations per tissue and type
valid_groups_tissue <- bind_rows(
    df_mean.dists_tissue %>% mutate(source = "mean"),
    df_median.dists_tissue %>% mutate(source = "median")
) %>%
    group_by(Tissue, Type, source) %>%
    summarise(n = n(), .groups = 'drop') %>%
    filter(n > 1) %>%
    split(.$source) %>%
    map(~select(.x, Tissue, Type))


# Perform t-tests for mean and median distances by tissue to compare between types (e.g., orthologs vs paralogs)
# Adjust p-values and apply the significance level function

# ------------------------------------------------------------------------

# Function to assign significance level based on p-value
significance_level <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")  
}

# ------------------------------------------------------------------------

# Perform t-tests for mean and median distances with error handling to catch any issues
# Perform both t-tests and Wilcoxon tests for tissue-specific comparisons
test_results_tissue <- tryCatch({
    # Median tests per tissue
    if (nrow(valid_groups_tissue$median) >= 2) {
        t_test_median_tissue <- df_median.dists_tissue %>%
            semi_join(valid_groups_tissue$median, by = c("Tissue", "Type")) %>%
            group_by(Tissue) %>%
            filter(!is.na(Distance)) %>%
            filter(n_distinct(Type) >= 2) %>%
            t_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(significance = sapply(p, significance_level),
                   analysis = "Median",
                   test_type = "t-test")
        
        wilcox_median_tissue <- df_median.dists_tissue %>%
            semi_join(valid_groups_tissue$median, by = c("Tissue", "Type")) %>%
            group_by(Tissue) %>%
            filter(!is.na(Distance)) %>%
            filter(n_distinct(Type) >= 2) %>%
            wilcox_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(significance = sapply(p, significance_level),
                   analysis = "Median",
                   test_type = "wilcox")
    }
    
    # Mean tests per tissue
    if (nrow(valid_groups_tissue$mean) >= 2) {
        t_test_mean_tissue <- df_mean.dists_tissue %>%
            semi_join(valid_groups_tissue$mean, by = c("Tissue", "Type")) %>%
            group_by(Tissue) %>%
            filter(!is.na(Distance)) %>%
            filter(n_distinct(Type) >= 2) %>%
            t_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(significance = sapply(p, significance_level),
                   analysis = "Mean",
                   test_type = "t-test")
        
        wilcox_mean_tissue <- df_mean.dists_tissue %>%
            semi_join(valid_groups_tissue$mean, by = c("Tissue", "Type")) %>%
            group_by(Tissue) %>%
            filter(!is.na(Distance)) %>%
            filter(n_distinct(Type) >= 2) %>%
            wilcox_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(significance = sapply(p, significance_level),
                   analysis = "Mean",
                   test_type = "wilcox")
    }
    
    # Combine all test results
    if (exists("t_test_median_tissue") && exists("t_test_mean_tissue") &&
        exists("wilcox_median_tissue") && exists("wilcox_mean_tissue")) {
        test_summary_tissue <- bind_rows(
            t_test_median_tissue, t_test_mean_tissue,
            wilcox_median_tissue, wilcox_mean_tissue
        )
        write.csv(test_summary_tissue, file.path(results_dir, "statistical_tests_summary_tissue.csv"), row.names = FALSE)
        message("Tissue-specific statistical tests summary exported to CSV")
    }
    
    list(
        t_test = list(
            median = if(exists("t_test_median_tissue")) t_test_median_tissue else NULL,
            mean = if(exists("t_test_mean_tissue")) t_test_mean_tissue else NULL
        ),
        wilcox = list(
            median = if(exists("wilcox_median_tissue")) wilcox_median_tissue else NULL,
            mean = if(exists("wilcox_mean_tissue")) wilcox_mean_tissue else NULL
        ),
        summary = if(exists("test_summary_tissue")) test_summary_tissue else NULL
    )
}, error = function(e) {
    message("Error in tissue-specific statistical tests: ", e$message)
    return(NULL)
})
