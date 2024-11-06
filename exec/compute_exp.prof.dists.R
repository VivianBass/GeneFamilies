require(GeneFamilies)
options(mc.cores = getMcCores())
library(dotenv)

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("USAGE: Rscript exec/compute_exp.prof.dists.R")

input.args <- commandArgs(trailingOnly = TRUE)
load(file.path(output_data_dir, "gene_families.RData"))                       
load(file.path(output_data_dir, "gene_groups_filtered.RData"))           
load(file.path(output_data_dir, "gene_expression.RData")) 

# would have to add the families also for distance calculation !?

# Function exp.prof.dists() sourced from:
source("R/compute_funks.R")

# Gene-Groups-Filtered:
# Initialize a character vector to store names of created objects

created_objects <- character()

# Calculate distances and add created objects to list if they exist
if (exists("con_orthologs_filtered.lst")) {
    con_orthologs.filtered.dists <- mclapply(con_orthologs_filtered.lst, exp.prof.dists)
    con_orthologs.filtered.dists.tissue <- mclapply(con_orthologs_filtered.lst, exp.prof.dists_tissue)
    created_objects <- c(created_objects, "con_orthologs.filtered.dists", "con_orthologs.filtered.dists.tissue")
}

if (exists("in_paralogs_filtered.lst")) {
    in_paralogs.filtered.dists <- mclapply(in_paralogs_filtered.lst, exp.prof.dists)
    in_paralogs.filtered.dists.tissue <- mclapply(in_paralogs_filtered.lst, exp.prof.dists_tissue)
    created_objects <- c(created_objects, "in_paralogs.filtered.dists", "in_paralogs.filtered.dists.tissue")
}

if (exists("out_paralogs_filtered.lst")) {
    out_paralogs.filtered.dists <- mclapply(out_paralogs_filtered.lst, exp.prof.dists)
    out_paralogs.filtered.dists.tissue <- mclapply(out_paralogs_filtered.lst, exp.prof.dists_tissue)
    created_objects <- c(created_objects, "out_paralogs.filtered.dists", "out_paralogs.filtered.dists.tissue")
}

if (exists("special_in_paralogs_filtered.lst")) {
    special_in_paralogs.filtered.dists <- mclapply(special_in_paralogs_filtered.lst, exp.prof.dists)
    special_in_paralogs.filtered.dists.tissue <- mclapply(special_in_paralogs_filtered.lst, exp.prof.dists_tissue)
    created_objects <- c(created_objects, "special_in_paralogs.filtered.dists", "special_in_paralogs.filtered.dists.tissue")
}

if (exists("special_out_paralogs_filtered.lst")) {
    special_out_paralogs.filtered.dists <- mclapply(special_out_paralogs_filtered.lst, exp.prof.dists)
    special_out_paralogs.filtered.dists.tissue <- mclapply(special_out_paralogs_filtered.lst, exp.prof.dists_tissue)
    created_objects <- c(created_objects, "special_out_paralogs.filtered.dists", "special_out_paralogs.filtered.dists.tissue")
}

# Save only the created objects
if (length(created_objects) > 0) {
    save(list = created_objects, file = file.path(output_data_dir, "exp.prof.dists_filtered.RData"))
    message("Objects have been saved to ", file.path(output_data_dir, "exp.prof.dists_filtered.RData"))
} else {
    message("No objects were created, so nothing was saved.")
}
