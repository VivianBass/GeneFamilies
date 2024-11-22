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

message("USAGE:  Rscript exec/load_expression_data.R <RPKM_counts_table.tsv>")

message("input.args[[1]]: <RPKM_counts_table.tsv>")
message("<RPKM_counts_table.tsv> expected to be TAB-Delimited")
message("<RPKM_counts_table.tsv> Header : \n", "id | tissue | expression")

# Parse input arguments
input.args <- commandArgs(trailingOnly = TRUE)

# read RPKM counts:
rpkm.rna.seq.counts <- read.table(input.args[[1]], sep = "\t", header = TRUE,
        check.names = FALSE, stringsAsFactors = FALSE) %>%
        select(id, tissue, expression) %>%
        mutate(expression = as.numeric(expression))

# create expression matrix, with all the tissues as Header 
expression_matrix <- rpkm.rna.seq.counts %>%
        pivot_wider(id_cols = id, names_from = tissue, values_from = expression,
        values_fn = list(expression = mean), values_fill = 0)

# normalize the expression matrix
rna.seq.exp.profils <- expression_matrix %>% rowwise() %>%
        mutate(row_sum = sum(c_across(-id), na.rm = TRUE)) %>%
        mutate(across(-c(id, row_sum), ~./row_sum)) %>%
        select(-row_sum) %>% ungroup()

# filter rna.seq.exp.profils for invalid or na values etc.
rna.seq.exp.profils <- rna.seq.exp.profils %>% rowwise() %>%
        filter(if_all(everything(), ~(!is.na(.) && . != "" && . != "NULL")))

# rename the expression profiles, and and map the protein identifiers if needed, to have matching with the names in the gene groups IDs
# basically add the Protein sequence ID to the expression profiles if required
# and !!! also add the species name to the expression profiles !!! extracted from the fasta files (larger_seq.fasta)
load("data/mapping_df.RData")
rna.seq.exp.profils <- rna.seq.exp.profils %>% rename(Parent_FBgn = id)
rna.seq.exp.profils <- rna.seq.exp.profils %>%
        left_join(mapping_df_unique, by = 'Parent_FBgn') %>%
        select(FBpp_ID, Parent_FBgn, Species, everything())
                                      
# Save results:
save(rna.seq.exp.profils, rpkm.rna.seq.counts, file = file.path(output_data_dir,"gene_expression.RData"))

write.table(rna.seq.exp.profils, file.path(output_data_dir, "RNA_Seq_RPKM_and_profiles.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

message("DONE")

