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
load(file.path(output_data_dir, "gene_expression.RData")) 

# Functions load_data_frame() & create_nested_list() sourced from:
source("R/load_data_funks.R")

# Load data frames using the load_data_frame() function from R/load_data_funks.R
con_orthologs <- load_data_frame(input.args[[1]])
in_paralogs <- load_data_frame(input.args[[2]])
out_paralogs <- load_data_frame(input.args[[3]])
special_in_paralogs <- load_data_frame(input.args[[4]])
special_out_paralogs <- load_data_frame(input.args[[5]])

# Save unfiltered data
save(
    con_orthologs, in_paralogs, out_paralogs, special_in_paralogs, special_out_paralogs,
    con_orthologs.lst, in_paralogs.lst, out_paralogs.lst, 
    special_in_paralogs.lst, special_out_paralogs.lst,
    file = file.path(output_data_dir, "gene_groups.RData")
)

# Filter out Data from the DataFrames which have no intersection with the provided expression data 
# this will reduce uneccessary computation time later on
# using filter_gene_pairs() from R/load_data_funks.R to reduce computation time.
con_orthologs_filtered <- filter_gene_pairs(con_orthologs, "Ortholog", rna.seq.exp.profils)
in_paralogs_filtered <- filter_gene_pairs(in_paralogs, "Paralog", rna.seq.exp.profils)
out_paralogs_filtered <- filter_gene_pairs(out_paralogs, "Paralog", rna.seq.exp.profils)
special_in_paralogs_filtered <- filter_gene_pairs(special_in_paralogs, "Paralog", rna.seq.exp.profils)
special_out_paralogs_filtered <- filter_gene_pairs(special_out_paralogs, "Paralog", rna.seq.exp.profils)

# Create nested lists with filtered data, using the create_nested_list() function from R/load_data_funks.R
con_orthologs_filtered.lst <- create_nested_list(con_orthologs_filtered, "Ortholog")
in_paralogs_filtered.lst <- create_nested_list(in_paralogs_filtered, "Paralog")
out_paralogs_filtered.lst <- create_nested_list(out_paralogs_filtered, "Paralog")
special_in_paralogs_filtered.lst <- create_nested_list(special_in_paralogs_filtered, "Paralog")
special_out_paralogs_filtered.lst <- create_nested_list(special_out_paralogs_filtered, "Paralog")


# Save filtered data
# Create list of objects to save
filtered_objects <- list()

# Check and add dataframes if they exist
if(exists("con_orthologs_filtered")) filtered_objects$con_orthologs_filtered <- con_orthologs_filtered
if(exists("in_paralogs_filtered")) filtered_objects$in_paralogs_filtered <- in_paralogs_filtered
if(exists("out_paralogs_filtered")) filtered_objects$out_paralogs_filtered <- out_paralogs_filtered
if(exists("special_in_paralogs_filtered")) filtered_objects$special_in_paralogs_filtered <- special_in_paralogs_filtered
if(exists("special_out_paralogs_filtered")) filtered_objects$special_out_paralogs_filtered <- special_out_paralogs_filtered

# Check and add nested lists if they exist
if(exists("con_orthologs_filtered.lst")) filtered_objects$con_orthologs_filtered.lst <- con_orthologs_filtered.lst
if(exists("in_paralogs_filtered.lst")) filtered_objects$in_paralogs_filtered.lst <- in_paralogs_filtered.lst
if(exists("out_paralogs_filtered.lst")) filtered_objects$out_paralogs_filtered.lst <- out_paralogs_filtered.lst
if(exists("special_in_paralogs_filtered.lst")) filtered_objects$special_in_paralogs_filtered.lst <- special_in_paralogs_filtered.lst
if(exists("special_out_paralogs_filtered.lst")) filtered_objects$special_out_paralogs_filtered.lst <- special_out_paralogs_filtered.lst

# Save only existing objects
do.call(save, c(names(filtered_objects), list(file = file.path(output_data_dir, "gene_groups_filtered.RData"))))

message("DONE")

