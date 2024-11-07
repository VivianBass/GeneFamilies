require(GeneFamilies)
options(mc.cores = getMcCores())
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(dotenv)

# Set up output directory, defined in the .env file 
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

cat("USAGE: Rscript exec/compute_exp.prof.dists.R")

# required data and files loaded from:
load(file.path(output_data_dir, "gene_expression.RData")) 
                    
# Functions sourced from:
source("R/compute_funks.R")

# would have to add the families also for distance calculation !?
# load(file.path(output_data_dir, "gene_families.RData"))


# load gene-groups datasets and catch the object names
load(file.path(output_data_dir, "gene_groups_filtered.RData")) 
loaded_objects <- ls()
gene_groups <- loaded_objects[grepl("_v\\.lst$", loaded_objects)]

# compute euclidean distances
for (group in gene_groups) {
    if (exists(group, envir = .GlobalEnv)) {
        data_object <- get(group)
        dist_name <- paste0(group, "_dists")
        assign(dist_name, mclapply(data_object, exp.prof.dists))
        created_objects <- c(created_objects, dist_name)
    }
}

# compute tissue-specific euclidean distances
for (group in gene_groups) {
    if (exists(group, envir = .GlobalEnv)) {
        data_object <- get(group)
        tissue_dist_name <- paste0(group, "_dists_tissue")
        assign(tissue_dist_name, mclapply(data_object, exp.prof.dists_tissue))
        created_objects <- c(created_objects, tissue_dist_name)
    }
}

save(list = created_objects, file = file.path(output_data_dir, "exp.prof.dists.RData"))
