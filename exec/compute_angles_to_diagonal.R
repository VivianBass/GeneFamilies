
library(dotenv)
# Define directories for output data and results using environment variables
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Create results directory if it doesn't exist
if(!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}

message("USAGE: Rscript exec/compute_exp.prof.dists_angles.R")

# Librarys for handling Dataframes, Lists etc. more efficiently
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(parallel)
                    
# functions sourced from:
source("R/angles_funks.R")
source("R/expression_funks.R")

# ------------------------------------------------------------------------

# load gene-groups datasets and catch the object names
load(file.path(output_data_dir, "gene_groups_filtered.RData")) 
loaded_objects <- ls()
gene_groups <- loaded_objects[grepl("_v\\.lst$", loaded_objects)]

# select your rna.seq.exp.profil data set and filter invalid Data, rows with NA etc
load(file.path(output_data_dir, "gene_expression.RData"))

tissues <- setdiff(colnames(rna.seq.exp.profils), c("FBpp_ID", "Species"))
# --------------------------------------------------------------------------------

angle_results <- list()

for(group in gene_groups) {

  gene_list <- get(group)
  df_name <- paste0(gsub("_v\\.lst$", "", group), ".expr.angle.diag.df")
  angle_results[[df_name]] <- calculate_angles(gene_list, rna.seq.exp.profils, tissues)
  assign(df_name, angle_results[[df_name]])
}

# Save the results
if(length(angle_results) == 0) {stop("No results to save")}

save(list = names(angle_results), file = file.path(output_data_dir, "exp.prof.dists_angles.RData"))






