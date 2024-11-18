require(GeneFamilies)
options(mc.cores = getMcCores())

cat("USAGE: Rscript exec/compute_exp.prof.dists.R")

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(dotenv)

# Set-up output directory, defined in the .env file 
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
                    
# functions sourced from:
source("R/compute_funks.R")

# load gene-groups datasets and catch the object names
load(file.path(output_data_dir, "gene_groups_filtered.RData")) 
loaded_objects <- ls()
gene_groups <- loaded_objects[grepl("_v\\.lst$", loaded_objects)]

# select your rna.seq.exp.profil data set and filter invalid Data, rows with NA etc
load(file.path(output_data_dir, "rna.seq.exp.profils_P_M_.RData"))

rna.seq.exp.profils <- rna.seq.exp.profils_P %>%
  distinct(Parent_FBgn, .keep_all = TRUE) %>%
  filter(!is.na(FBpp_ID) & FBpp_ID != "NA" & FBpp_ID != "")

tissues <- setdiff(colnames(rna.seq.exp.profils), c("FBpp_ID", "Parent_FBgn", "Species"))

rna.seq.exp.profils <- rna.seq.exp.profils %>%
  filter(rowSums(across(all_of(tissues), 
  ~(. == "Invalid Number" | is.na(.) | . == "NA" | . == "" | . == "NaN" |
  . == "missing"))) != length(tissues))

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
cat("DONE")
