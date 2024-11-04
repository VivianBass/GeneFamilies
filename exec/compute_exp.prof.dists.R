require(GeneFamilies)
options(mc.cores = getMcCores())
library(dotenv)

library(dplyr)
library(tidyr)
library(tibble)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("USAGE: Rscript exec/compute_exp.prof.dists.R")

input.args <- commandArgs(trailingOnly = TRUE)
load(file.path(output_data_dir, "gene_families.RData"))                       
load(file.path(output_data_dir, "gene_groups.RData"))              
load(file.path(output_data_dir, "gene_expression.RData"))     

# Function exp.prof.dists() sourced from:
source("R/compute_funks.R")

# Gene-Groups:
con_orthologs.dists <- mclapply(con_orthologs.lst, exp.prof.dists)
con_orthologs.dists.tissue <- mclapply(con_orthologs.lst, exp.prof.dists, per.tissue = TRUE)

in_paralogs.dists <- mclapply(in_paralogs.lst, exp.prof.dists)
in_paralogs.dists.tissue <- mclapply(in_paralogs.lst, exp.prof.dists, per.tissue = TRUE)

out_paralogs.dists <- mclapply(out_paralogs.lst, exp.prof.dists)
out_paralogs.dists.tissue <- mclapply(out_paralogs.lst, exp.prof.dists, per.tissue = TRUE)

special_in_paralogs.dists <- mclapply(special_in_paralogs.lst, exp.prof.dists)
special_in_paralogs.dists.tissue <- mclapply(special_in_paralogs.lst, exp.prof.dists, per.tissue = TRUE)

special_out_paralogs.dists <- mclapply(special_out_paralogs.lst, exp.prof.dists)
special_out_paralogs.dists.tissue <- mclapply(special_out_paralogs.lst, exp.prof.dists, per.tissue = TRUE)

# Gene-Families:
non.singleton.fams <- families.df$id[which(families.df$size > 1)]
families.exp.prof.dists <- mclapply(families.lst[non.singleton.fams], exp.prof.dists_2)
families.exp.prof.dists.tissue <- mclapply(families.lst[non.singleton.fams], 
                                        exp.prof.dists_2, per.tissue = TRUE)

save(con_orthologs.dists, con_orthologs.dists.tissue,
     in_paralogs.dists, in_paralogs.dists.tissue,
     out_paralogs.dists, out_paralogs.dists.tissue,
     special_in_paralogs.dists, special_in_paralogs.dists.tissue,
     special_out_paralogs.dists, special_out_paralogs.dists.tissue,
     families.exp.prof.dists, families.exp.prof.dists.tissue,
     file = file.path(output_data_dir, "exp.prof.dists.RData"))


message("DONE")

