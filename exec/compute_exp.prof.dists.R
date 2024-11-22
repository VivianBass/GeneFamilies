require(GeneFamilies)
options(mc.cores = getMcCores())
library(parallel)

library(dotenv)
# Define directories for output data and results using environment variables
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Librarys for handling Dataframes, Lists etc. more efficiently
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

message("USAGE: Rscript exec/compute_exp.prof.dists.R")
             
# functions sourced from:
source("R/compute_funks.R")

# load gene-groups datasets and catch the object names
load(file.path(output_data_dir, "gene_groups_filtered.RData")) 
loaded_objects <- ls()
gene_groups <- loaded_objects[grepl("_v\\.lst$", loaded_objects)]

# select your rna.seq.exp.profil data set and filter invalid Data, rows with NA etc
load(file.path(output_data_dir, "gene_expression.RData"))

rna.seq.exp.profils <- rna.seq.exp.profils %>%
    distinct(Parent_FBgn, .keep_all = TRUE) %>%
    filter(!is.na(FBpp_ID) & FBpp_ID != "NA" & FBpp_ID != "")

tissues <- setdiff(colnames(rna.seq.exp.profils), c("FBpp_ID", "Parent_FBgn", "Species"))

# filter rna.seq.exp.profils for invalid or na values etc
rna.seq.exp.profils <- rna.seq.exp.profils %>%
    filter(rowSums(across(all_of(tissues), 
    ~(. == "Invalid Number" | is.na(.) | . == "NA" | . == "" | . == "NaN" |
    . == "missing"))) != length(tissues))

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
message("DONE")



