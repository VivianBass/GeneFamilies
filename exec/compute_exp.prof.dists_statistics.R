require(GeneFamilies)
options(mc.cores = getMcCores())
library(dotenv)

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("USAGE: Rscript exec/compute_exp.prof.dists_statistics.R")

load(file.path(output_data_dir, "exp.prof.dists_unfiltered.RData"))
load(file.path(output_data_dir, "exp.prof.dists_families.RData"))
load(file.path(output_data_dir, "exp.prof.dists_filtered.RData"))

# Functions sourced from: source("R/compute_funks.R")
# add_to_list_if_exists()
# calculate_exp.prof.dists.statistics() & calculate_exp.prof.dists.tissue.statistics() 

source("R/compute_funks.R")

# Set up log file
log_file <- file.path("exp.prof.dists_statistics_log.txt")
sink(log_file, append = TRUE)
cat("Starting statistics computation...\n\n")

# Calculate statistics for filtered data

filtered_stats <- list()
filtered_stats <- add_to_list_if_exists(filtered_stats, "df_median_mean_con_orthologs", "con_orthologs.filtered.dists")
filtered_stats <- add_to_list_if_exists(filtered_stats, "df_median_mean_in_paralogs", "in_paralogs.filtered.dists")
filtered_stats <- add_to_list_if_exists(filtered_stats, "df_median_mean_out_paralogs", "out_paralogs.filtered.dists")
filtered_stats <- add_to_list_if_exists(filtered_stats, "df_median_mean_special_in_paralogs", "special_in_paralogs.filtered.dists")
filtered_stats <- add_to_list_if_exists(filtered_stats, "df_median_mean_special_out_paralogs", "special_out_paralogs.filtered.dists")

# Save filtered statistics and log
if (length(filtered_stats) > 0) {
    save(filtered_stats, file = file.path(output_data_dir, "exp.prof.dists_filtered_statistics.RData"))
    cat("Filtered statistics saved successfully.\n\n")
} else {
    cat("No filtered statistics objects to save.\n\n")
}

cat("Statistics computation complete.\n")
sink()  # Close the sink

message(" DONE. Check the log file for details: ", log_file)


