require(GeneFamilies)
options(mc.cores = getMcCores())
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(dotenv)

# Set up output directory, defined in the .env file 
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("USAGE: Rscript exec/load_gene_groups_data.R input.args ...")
message("PURPOSE: This R script loads five different gene groups as stated in the messages below") 
message("All files should follow the same header naming convention as stated in Header-Type messages below. ") 

message("input.args[[1]]:  path/2/<conserved_orthologs.tsv>")
message("input.args[[2]]:  path/2/<in_paralogs.tsv>")
message("input.args[[3]]:  path/2/<out_paralogs.tsv>")
message("input.args[[4]]:  path/2/<special_in_paralogs.tsv>")
message("input.args[[5]]:  path/2/<special_out_paralogs.tsv>")

input.args <- commandArgs(trailingOnly = TRUE)

# required data and files loaded from:
load(file.path(output_data_dir, "gene_expression_flybase.RData")) 

# Functions load_data_frame() & create_nested_list() sourced from:
source("R/load_data_funks.R")

# define input file names and their corresponding types
input_files <- list(
    con_orthologs = list(path = input.args[[1]], type = "Ortholog"),
    in_paralogs = list(path = input.args[[2]], type = "Paralog"),
    out_paralogs = list(path = input.args[[3]], type = "Paralog"),
    special_in_paralogs = list(path = input.args[[4]], type = "Paralog"),
    special_out_paralogs = list(path = input.args[[5]], type = "Paralog")
)

# track created objects
created_objects <- character()

# load data frames and create nested lists
for (name in names(input_files)) {
    # load data frame
    df_name <- name
    assign(df_name, load_data_frame(input_files[[name]]$path))
    created_objects <- c(created_objects, df_name)
    
    # create nested list
    lst_name <- paste0(name, ".lst")
    assign(lst_name, create_nested_list(get(df_name), input_files[[name]]$type))
    created_objects <- c(created_objects, lst_name)
}

save(list = created_objects, file = file.path(output_data_dir, "gene_groups_unfiltered.RData"))


# Filtering process: using either filter_v1() or filter_v2() functions to filter out genes
# that are not present in rna.seq.exp.profiles and therefore lack expression values.
# Removing these genes helps reduce congestion in subsequent computations.
# Both filters use the intersect() method to retain only relevant genes.
# For greater accuracy, use filter_v2() to filter by two columns

filtered_objects <- list()

for (name in names(input_files)) {
    tryCatch({
        filtered_name <- paste0(name, "_v")
        if (exists(name, envir = .GlobalEnv)) {
            assign(filtered_name, filter_v1(get(name), rna.seq.exp.profils))
            filtered_objects[[filtered_name]] <- get(filtered_name)
            cat(sprintf("Object '%s' added successfully.\n", filtered_name))
            
            # Only create nested list if filtering was successful
            if (nrow(get(filtered_name)) > 0) {
                filtered_lst_name <- paste0(filtered_name, ".lst")
                assign(filtered_lst_name, create_nested_list(get(filtered_name), input_files[[name]]$type))
                filtered_objects[[filtered_lst_name]] <- get(filtered_lst_name)
                cat(sprintf("Object '%s' added successfully.\n", filtered_lst_name))
            }
        }
    }, error = function(e) {cat(sprintf("Error processing %s: %s\n", name, e$message))})
}

# Save filtered objects
if (length(names(filtered_objects)) > 0) {
    save(list = names(filtered_objects), 
         file = file.path(output_data_dir, "gene_groups_filtered.RData"))
}
message("DONE")



