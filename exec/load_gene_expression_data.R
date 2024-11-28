
library(dotenv)
# Define directories for output data and results using environment variables
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Create results directory if it doesn't exist
if(!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}

message("USAGE:  Rscript exec/load_expression_data.R <RPKM_counts_table.tsv>")

# Parse input arguments
input.args <- commandArgs(trailingOnly = TRUE)

message("input.args[[1]]: <RPKM_counts_table.tsv>")
message("<RPKM_counts_table.tsv> expected to be TAB-Delimited")
message("<RPKM_counts_table.tsv> Header : \n", "id | tissue | expression")

# Librarys for handling Dataframes, Lists etc. more efficiently
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(parallel)

# ------------------------------------------------------------------------

### Log2 Transformation Effects:
### - Values between 0 and 1 become negative (closer to 0 -> closer to 0, farther -> large negatives)
### - Values between 1 and 2 transform to small positives (closer to 2 -> closer to 1)
### - Values greater than 2 grow into larger positives (e.g., 2 -> 1, 4 -> 2, 10 -> ~3.32)

# read RPKM counts:
rpkm.rna.seq.counts <- read.table(input.args[[1]], sep = "\t", header = TRUE,
        check.names = FALSE, stringsAsFactors = FALSE) %>%
        select(id, tissue, expression) %>%
        mutate(expression = as.numeric(expression))

rpkm.rna.seq.counts_log2 <- rpkm.rna.seq.counts %>%
        mutate(expression = log2(expression + 1))

# create expression matrix, with all the tissues as Header 
expression_matrix <- rpkm.rna.seq.counts %>%
        pivot_wider(id_cols = id, names_from = tissue, values_from = expression,
        values_fn = list(expression = mean), values_fill = 0)

expression_matrix_log2 <- rpkm.rna.seq.counts_log2 %>%
        pivot_wider(id_cols = id, names_from = tissue, values_from = expression,
        values_fn = list(expression = mean), values_fill = 0)

# normalize the expression matrix
rna.seq.exp.profils <- expression_matrix %>% rowwise() %>%
        mutate(row_sum = sum(c_across(-id), na.rm = TRUE)) %>%
        mutate(across(-c(id, row_sum), ~./row_sum)) %>%
        select(-row_sum) %>% ungroup()

rna.seq.exp.profils_log2 <- expression_matrix_log2 %>% rowwise() %>%
        mutate(row_sum = sum(c_across(-id), na.rm = TRUE)) %>%
        mutate(across(-c(id, row_sum), ~./row_sum)) %>%
        select(-row_sum) %>% ungroup()

# filter rna.seq.exp.profils for invalid or na values etc.
# and rename the first column to FBpp_ID (column containing Gene IDs)
rna.seq.exp.profils <- rna.seq.exp.profils %>%
        rename(FBpp_ID = 1) %>%
        rowwise() %>%
        filter(if_all(everything(), ~(!is.na(.) && . != "" && . != "NULL"))) %>%
        distinct(FBpp_ID, .keep_all = TRUE)

rna.seq.exp.profils_log2 <- rna.seq.exp.profils_log2 %>%
        rename(FBpp_ID = 1) %>%
        rowwise() %>%
        filter(if_all(everything(), ~(!is.na(.) && . != "" && . != "NULL"))) %>%
        distinct(FBpp_ID, .keep_all = TRUE)

# filter rna.seq.exp.profils for invalid or na values etc
tissues <- setdiff(colnames(rna.seq.exp.profils), c("FBpp_ID", "Species"))

rna.seq.exp.profils <- rna.seq.exp.profils %>%
        filter(rowSums(across(all_of(tissues), 
        ~(. == "Invalid Number" | is.na(.) | . == "NA" | . == "" | . == "NaN" |
        . == "missing"))) != length(tissues))

tissues_log2 <- setdiff(colnames(rna.seq.exp.profils_log2), c("FBpp_ID", "Species"))

rna.seq.exp.profils_log2 <- rna.seq.exp.profils_log2 %>%
        filter(rowSums(across(all_of(tissues_log2), 
        ~(. == "Invalid Number" | is.na(.) | . == "NA" | . == "" | . == "NaN" |
        . == "missing"))) != length(tissues))

# Merge expression profiles with mapping data, remove Parent_FBgn column,
# and reorder columns to have FBpp_ID and Species first
load(file.path(output_data_dir, "mapping_df.RData"))

rna.seq.exp.profils <- rna.seq.exp.profils %>%
        left_join(mapping_df, by = "FBpp_ID") %>%
        select(-Parent_FBgn) %>% 
        select(FBpp_ID, Species, everything())

rna.seq.exp.profils_log2 <- rna.seq.exp.profils_log2 %>%
        left_join(mapping_df, by = "FBpp_ID") %>%
        select(-Parent_FBgn) %>% 
        select(FBpp_ID, Species, everything())
                         
# Save results:
save(rna.seq.exp.profils, rpkm.rna.seq.counts,
     rna.seq.exp.profils_log2, rpkm.rna.seq.counts_log2, 
     file = file.path(output_data_dir,"gene_expression.RData"))

write.table(rna.seq.exp.profils, file.path(output_data_dir, "RNA_Seq_RPKM_and_profiles.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(rna.seq.exp.profils_log2, file.path(output_data_dir, "RNA_Seq_RPKM_and_profiles_log2.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)


