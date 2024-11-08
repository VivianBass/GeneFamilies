require(GeneFamilies)
options(mc.cores = getMcCores())
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(dotenv)

# Set-up output directory, defined in the .env file 
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

cat("USAGE: Rscript exec/compute_exp.prof.dists_statistics.R")

# functions sourced from:
source("R/compute_funks.R")

# Automatically sort the loaded gene groups data into regular and tissue datasets
load(file.path(output_data_dir, "exp.prof.dists.RData"))
loaded_objects <- ls()
# This function validates loaded data, and returns a vector of valid names.
valid_data_names <- validate_data(loaded_objects, "(.lst_dists$)")
valid_data_names_tissue <- validate_data(loaded_objects, "(.lst_dists_tissue$)")

# Track created dataframe names
created_dfs <- character()

# compute mean/median statistics
for (name in valid_data_names) {
    if (exists(name, envir = .GlobalEnv)) {
        data_object <- get(name)
        df_name <- paste0(name, "_stats")
        assign(df_name, calculate_exp.prof.dists.statistics(data_object))
        created_dfs <- c(created_dfs, df_name)
    }
}

# compute mean/median tissue-specific statistics
for (name in valid_data_names_tissue ) {
    if (exists(name, envir = .GlobalEnv)) {
        data_object <- get(name)
        df_name <- paste0(name, "_stats")
        assign(df_name, calculate_exp.prof.dists.tissue.statistics(data_object))
        created_dfs <- c(created_dfs, df_name)
    }
}

save(list = created_dfs, file = file.path(output_data_dir, "exp.prof.dists_statistics.RData"))

