
require(GeneFamilies)
options(mc.cores = getMcCores())

cat("USAGE: Rscript exec/compute_exp.prof.dists_angles.R")

library(parallel)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(dotenv)

# Set-up output directory, defined in the .env file 
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
                    
# functions sourced from:
source("R/expression_funks.R")

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

# calculate exp.prof.dists angles to diagonal for con_orthologs 
con_orthologs.genes <- unlist(con_orthologs_v.lst)
con_orthologs.expr <- intersect(con_orthologs.genes, rna.seq.exp.profils$FBpp_ID)

con_orthologs.expr.angle.diag.df <- data.frame(FBpp_ID = orths.expr, angle.diag = as.numeric(mclapply(orths.expr, 
    function(x) {
        cosDiag(rna.seq.exp.profils[which(rna.seq.exp.profils$FBpp_ID == 
            x), tissues])/sqrt(2)
    })), stringsAsFactors = FALSE)

con_orthologs.expr.angle.diag.df <- con_orthologs.expr.angle.diag.df %>%
    filter(!is.na(angle.diag) & angle.diag != "NA" & angle.diag != "")

# calculate exp.prof.dists angles to diagonal for in_paralogs.lst 
in_paralogs.genes <- unlist(in_paralogs_v.lst)
in_paralogs.expr <- intersect(in_paralogs.genes, rna.seq.exp.profils$FBpp_ID)

in_paralogs.expr.angle.diag.df <- data.frame(FBpp_ID = in_paralogs.expr, angle.diag = as.numeric(mclapply(in_paralogs.expr, 
    function(x) {
        cosDiag(rna.seq.exp.profils[which(rna.seq.exp.profils$FBpp_ID == 
            x), tissues])/sqrt(2)
    })), stringsAsFactors = FALSE)

in_paralogs.expr.angle.diag.df <- in_paralogs.expr.angle.diag.df %>%
    filter(!is.na(angle.diag) & angle.diag != "NA" & angle.diag != "")

# calculate exp.prof.dists angles to diagonal for out_paralogs.lst
out_paralogs.genes <- unlist(out_paralogs_v.lst)
out_paralogs.expr <- intersect(out_paralogs.genes, rna.seq.exp.profils$FBpp_ID)

out_paralogs.expr.angle.diag.df <- data.frame(FBpp_ID = out_paralogs.expr, angle.diag = as.numeric(mclapply(out_paralogs.expr, 
    function(x) {
        cosDiag(rna.seq.exp.profils[which(rna.seq.exp.profils$FBpp_ID == 
            x), tissues])/sqrt(2)
    })), stringsAsFactors = FALSE)

out_paralogs.expr.angle.diag.df <- out_paralogs.expr.angle.diag.df %>%
    filter(!is.na(angle.diag) & angle.diag != "NA" & angle.diag != "")

# calculate exp.prof.dists angles to diagonal for special_in_paralogs.lst 
special_in_paralogs.genes <- unlist(special_in_paralogs_v.lst)
special_in_paralogs.expr <- intersect(special_in_paralogs.genes, rna.seq.exp.profils$FBpp_ID)

special_in_paralogs.expr.angle.diag.df <- data.frame(FBpp_ID = special_in_paralogs.expr, angle.diag = as.numeric(mclapply(special_in_paralogs.expr, 
    function(x) {
        cosDiag(rna.seq.exp.profils[which(rna.seq.exp.profils$FBpp_ID == 
            x), tissues])/sqrt(2)
    })), stringsAsFactors = FALSE)

special_in_paralogs.expr.angle.diag.df <- special_in_paralogs.expr.angle.diag.df %>%
    filter(!is.na(angle.diag) & angle.diag != "NA" & angle.diag != "")

# calculate exp.prof.dists angles to diagonal for special_out_paralogs.lst 
special_out_paralogs.genes <- unlist(special_out_paralogs_v.lst)
special_out_paralogs.expr <- intersect(special_out_paralogs.genes, rna.seq.exp.profils$FBpp_ID)

special_out_paralogs.expr.angle.diag.df <- data.frame(FBpp_ID = special_out_paralogs.expr, angle.diag = as.numeric(mclapply(special_out_paralogs.expr, 
    function(x) {
        cosDiag(rna.seq.exp.profils[which(rna.seq.exp.profils$FBpp_ID == 
            x), tissues])/sqrt(2)
    })), stringsAsFactors = FALSE)

special_out_paralogs.expr.angle.diag.df <- special_out_paralogs.expr.angle.diag.df %>%
    filter(!is.na(angle.diag) & angle.diag != "NA" & angle.diag != "")

# --------------------------------------------------------------------------------

# Function to calculate angles for a gene group
calculate_angles <- function(genes, rna.seq.exp.profils, tissues) {
  genes.expr <- intersect(unlist(genes), rna.seq.exp.profils$FBpp_ID)
  
  if(length(genes.expr) == 0) {
    warning(paste("No matching genes found for group:", group))
    return(data.frame())
  }

  expr.angle.diag.df <- data.frame(
    FBpp_ID = genes.expr,
    angle.diag = as.numeric(mclapply(genes.expr, function(x) {
      cosDiag(rna.seq.exp.profils[which(rna.seq.exp.profils$FBpp_ID == x), tissues])/sqrt(2)
    })),
    stringsAsFactors = FALSE
  )

  expr.angle.diag.df %>%
    filter(!is.na(angle.diag) & angle.diag != "NA" & angle.diag != "")
}

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

