require(GeneFamilies)
options(mc.cores = getMcCores())
library(dotenv)

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

# Set up output directory, defined in the .env file 
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

# Display usage and purpose information
message("USAGE: Rscript exec/load_gene_groups_data.R input.args ...")
message("PURPOSE: This R script loads five different gene groups of pairwise orthologs and paralogs data, including one conserved ortholog file and four different paralog files. All files should follow the same header naming convention as stated in Header-Type messages below.")

message("Header-Type-Orthologs: Family | Gene | Gene_species | Ortholog | Ortholog_species")
message("input.args[[1]]:  path/2/<conserved_orthologs.tsv>")
message("Header-Type-Paralogs:  Family | Gene | Gene_species | Paralog  | Paralog_species")
message("input.args[[2]]:  path/2/<in_paralogs.tsv>")
message("input.args[[3]]:  path/2/<out_paralogs.tsv>")
message("input.args[[4]]:  path/2/<special_in_paralogs.tsv>")
message("input.args[[5]]:  path/2/<special_out_paralogs.tsv>")

# Parse input arguments
input.args <- commandArgs(trailingOnly = TRUE)

# Functions load_data_frame() & create_nested_list() sourced from:
source("R/loading_section_funks.R")

# Load data frames using the load_data_frame() function from R/loading_section_funks.R
con_orthologs <- load_data_frame(input.args[[1]])
in_paralogs <- load_data_frame(input.args[[2]])
out_paralogs <- load_data_frame(input.args[[3]])
special_in_paralogs <- load_data_frame(input.args[[4]])
special_out_paralogs <- load_data_frame(input.args[[5]])

# Create nested lists using the create_nested_list() function from R/loading_section_funks.R
con_orthologs.lst <- create_nested_list(con_orthologs, "Ortholog")
in_paralogs.lst <- create_nested_list(in_paralogs, "Paralog")
out_paralogs.lst <- create_nested_list(out_paralogs, "Paralog")
special_in_paralogs.lst <- create_nested_list(special_in_paralogs, "Paralog")
special_out_paralogs.lst <- create_nested_list(special_out_paralogs, "Paralog")


save(con_orthologs, con_orthologs.lst, in_paralogs, in_paralogs.lst,
     out_paralogs, out_paralogs.lst, special_in_paralogs, special_in_paralogs.lst,
     special_out_paralogs, special_out_paralogs.lst, 
     file = file.path(output_data_dir, "gene_groups.RData"))

message("DONE")