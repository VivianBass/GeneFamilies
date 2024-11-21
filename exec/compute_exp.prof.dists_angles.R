require(GeneFamilies)
options(mc.cores = getMcCores())
library(parallel)

message("USAGE: Rscript exec/compute_exp.prof.dists_angles.R")

library(dotenv)
# Define directories for output data and results using environment variables
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Librarys for handling Dataframes, Lists etc. more efficiently
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
                    
# functions sourced from:
source("R/angles_funks.R")
source("R/expression_funks.R")

# ------------------------------------------------------------------------

# load gene-groups datasets and catch the object names
load(file.path(output_data_dir, "gene_groups_filtered.RData")) 
loaded_objects <- ls()
gene_groups <- loaded_objects[grepl("_v\\.lst$", loaded_objects)]

# select your rna.seq.exp.profil data set and filter invalid Data, rows with NA etc
load(file.path(output_data_dir, "rna.seq.exp.profils_P_M_.RData"))

rna.seq.exp.profils <- rna.seq.exp.profils_M %>%
    distinct(Parent_FBgn, .keep_all = TRUE) %>%
    filter(!is.na(FBpp_ID) & FBpp_ID != "NA" & FBpp_ID != "")

tissues <- setdiff(colnames(rna.seq.exp.profils), c("FBpp_ID", "Parent_FBgn", "Species"))

# --------------------------------------------------------------------------------

angle_results <- list()

# Loop through each gene group
for(group in gene_groups) {
  gene_list <- get(group)
  
  # Create dataframe name dynamically 
  df_name <- paste0(gsub("_v\\.lst$", "", group), ".expr.angle.diag.df")
  
  # Store result in angle_results list
  angle_results[[df_name]] <- calculate_angles(gene_list, rna.seq.exp.profils, tissues)
  
  # Also assign to global environment
  assign(df_name, angle_results[[df_name]])
}

# Save the results
if(length(angle_results) == 0) {
  stop("No results to save")
}

save(list = names(angle_results), file = file.path(output_data_dir, "exp.prof.dists_angles.RData"))
message("DONE")





