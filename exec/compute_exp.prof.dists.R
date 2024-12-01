
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
created_objects_euclidean <- c()

# compute euclidean distances
for (group in gene_groups) {
    if (exists(group, envir = .GlobalEnv)) {
        data_object <- get(group)
        dist_name <- paste0(group, "_dists")
        assign(dist_name, mclapply(data_object, exp.prof.dists))
        created_objects_euclidean <- c(created_objects_euclidean, dist_name)
    }
}

# compute tissue-specific euclidean distances
for (group in gene_groups) {
    if (exists(group, envir = .GlobalEnv)) {
        data_object <- get(group)
        tissue_dist_name <- paste0(group, "_dists_tissue")
        assign(tissue_dist_name, mclapply(data_object, exp.prof.dists_tissue))
        created_objects_euclidean <- c(created_objects_euclidean, tissue_dist_name)
    }
}

save(list = created_objects_euclidean, file = file.path(output_data_dir, "exp.prof.dists.RData"))

# ------------------------------------------------------------------------
# using log2 transformed expression values for distance calculation

# Compute Euclidean distances with log2 transformed values
created_objects_euclidean_log2 <- c()

for (group in gene_groups) {
    if (exists(group, envir = .GlobalEnv)) {
        data_object <- get(group)
        dist_name <- paste0(group, "_dists_log2")
        
        # Pass the log2-transformed data frame to the function
        assign(dist_name, mclapply(data_object, exp.prof.dists, expression.profiles = rna.seq.exp.profils_log2))
        
        created_objects_euclidean_log2 <- c(created_objects_euclidean_log2, dist_name)
    }
}

# Compute tissue-specific Euclidean distances with log2 transformed values
for (group in gene_groups) {
    if (exists(group, envir = .GlobalEnv)) {
        data_object <- get(group)
        tissue_dist_name <- paste0(group, "_dists_tissue_log2")
        
        # Pass the log2-transformed data frame to the function
        assign(tissue_dist_name, mclapply(data_object, exp.prof.dists_tissue, expression.profiles = rna.seq.exp.profils_log2))
        
        created_objects_euclidean_log2 <- c(created_objects_euclidean_log2, tissue_dist_name)
    }
}

# Save the log2-transformed distance objects
save(list = created_objects_euclidean_log2, file = file.path(output_data_dir, "exp.prof.dists.log2.RData"))

# ------------------------------------------------------------------------
# using cosine angles for distance calculation

# Initialize created_objects vector before the loops
created_objects_angles <- c()

# compute cosine angles distances
for (group in gene_groups) {
    if (exists(group, envir = .GlobalEnv)) {
        data_object <- get(group)
        dist_name <- paste0(group, "_cos_angles_dists")
        assign(dist_name, mclapply(data_object, exp.prof.angles))
        created_objects_angles <- c(created_objects_angles, dist_name)
    }
}

# compute cosine angles distances
#for (group in gene_groups) {
#    if (exists(group, envir = .GlobalEnv)) {
#        data_object <- get(group)
#        dist_name <- paste0(group, "_cos_angles_dists")
#        assign(dist_name, mclapply(data_object, exp.prof.dists, dist.method = "angle"))
#        created_objects_angles <- c(created_objects_angles, dist_name)
#    }
#}

save(list = created_objects_angles , file = file.path(output_data_dir, "exp.prof.angles.RData"))

# ------------------------------------------------------------------------
# using cosine angles for distance calculation with log2 transformed values

created_objects_angles_log2 <- c()

# compute cosine angles distances
for (group in gene_groups) {
    if (exists(group, envir = .GlobalEnv)) {
        data_object <- get(group)
        dist_name <- paste0(group, "_cos_angles_dists_log2")
        assign(dist_name, mclapply(data_object, exp.prof.angles, expression.profiles = rna.seq.exp.profils_log2))
        created_objects_angles_log2 <- c(created_objects_angles_log2, dist_name)
    }
}

#for (group in gene_groups) {
#    if (exists(group, envir = .GlobalEnv)) {
#        data_object <- get(group)
#        dist_name <- paste0(group, "_cos_angles_dists_log2")
#        assign(dist_name, mclapply(data_object, exp.prof.dists, expression.profiles = rna.seq.exp.profils_log2, dist.method = "angle"))
#        created_objects_angles <- c(created_objects_angles, dist_name)
#    }
#}

save(list = created_objects_angles_log2, file = file.path(output_data_dir, "exp.prof.angles.log2.RData"))
