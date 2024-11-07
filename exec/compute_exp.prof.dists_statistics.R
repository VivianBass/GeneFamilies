require(GeneFamilies)
options(mc.cores = getMcCores())
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(dotenv)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

cat("USAGE: Rscript exec/compute_exp.prof.dists_statistics.R\n")

# Load functions:
source("R/compute_funks.R")

# would have to add the families also for distance calculation !?
# load(file.path(output_data_dir, "gene_families.RData"))

# Automatically sort the loaded gene groups data into regular and tissue datasets
load(file.path(output_data_dir, "exp.prof.dists.RData"))
loaded_objects <- ls()
data_names <- loaded_objects[grepl("\\.filtered\\.dists$", loaded_objects)]
data_names_tissue <- loaded_objects[grepl("\\.filtered\\.dists\\.tissue$", loaded_objects)]

# Track created dataframe names
created_dfs <- character()

# compute mean/median statistics
for (name in data_names) {
    if (exists(name, envir = .GlobalEnv)) {
        data_object <- get(name)
        df_name <- paste0(name, "_stats")
        assign(df_name, calculate_exp.prof.dists.statistics(data_object))
        created_dfs <- c(created_dfs, df_name)
    }
}

# compute mean/median tissue-specific statistics
for (name in data_names_tissue) {
    if (exists(name, envir = .GlobalEnv)) {
        data_object <- get(name)
        df_name <- paste0(name, "_stats")
        assign(df_name, calculate_exp.prof.dists.tissue.statistics(data_object))
        created_dfs <- c(created_dfs, df_name)
    }
}

save(list = created_dfs, file = file.path(output_data_dir, "exp.prof.dists_statistics.RData"))

