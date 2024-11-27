
library(dotenv)
# Define directories for output data and results using environment variables
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Create results directory if it doesn't exist
if(!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}

message("USAGE: Rscript exec/compute_exp.prof.dists.R")

# Librarys for handling Dataframes, Lists etc. more efficiently
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(parallel)
             
# functions sourced from:
source("R/compute_funks.R")

# ------------------------------------------------------------------------

# load gene-groups datasets and catch the object names
load(file.path(output_data_dir, "gene_groups_filtered.RData")) 
loaded_objects <- ls()
gene_groups <- loaded_objects[grepl("_v\\.lst$", loaded_objects)]

# select your rna.seq.exp.profil data set and filter invalid Data, rows with NA etc
load(file.path(output_data_dir, "gene_expression.RData"))

# ------------------------------------------------------------------------

# Initialize created_objects vector before the loops
created_objects <- c()

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

# ------------------------------------------------------------------------


# other distance methods (angles)

# created_objects <- c()
# compute euclidean distances
# for (group in gene_groups) {
#    if (exists(group, envir = .GlobalEnv)) {
#        data_object <- get(group)
#        dist_name <- paste0(group, "_log2_dists")
#        assign(dist_name, mclapply(data_object, exp.prof.dists_log2 ))
#        created_objects <- c(created_objects, dist_name)
#    }
#}

#created_objects <- c()
# compute euclidean distances
#for (group in gene_groups) {
#    if (exists(group, envir = .GlobalEnv)) {
#        data_object <- get(group)
#        dist_name <- paste0(group, "_cosine_dists")
#        assign(dist_name, mclapply(data_object, exp.prof_cosine ))
#        created_objects <- c(created_objects, dist_name)
#    }
#}



