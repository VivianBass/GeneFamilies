require(GeneFamilies)
options(mc.cores = getMcCores())
library(dotenv)
library(parallel)

library(dplyr)
library(tidyr)
library(tibble)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("USAGE:  Rscript exec/load_expression_data.R <RPKM_counts_table.tsv>")

message("input.args[[1]]: <RPKM_counts_table.tsv>")
message("<RPKM_counts_table.tsv> expected to be TAB-Delimited")
message("<RPKM_counts_table.tsv> Header : \n", "id / tissue / expression")

input.args <- commandArgs(trailingOnly = TRUE)

# read RPKM counts:
rpkm.rna.seq.counts <- read.table(input.args[[1]], 
        sep = "\t", header = TRUE, fill = TRUE,
        stringsAsFactors = FALSE, comment.char = "", quote = "", na.strings = "") %>%
        select(id, tissue, expression) %>%
        mutate(expression = as.numeric(expression))

# create expression matrix, with all the tissues as Header 
expression_matrix <- rpkm.rna.seq.counts %>%
        pivot_wider(id_cols = id, names_from = tissue, values_from = expression,
        values_fn = list(expression = mean), values_fill = 0)

expression_matrix <- expression_matrix %>% select(-`NA`)

# normalize the expression matrix
rna.seq.exp.profils <- expression_matrix %>% rowwise() %>%
        mutate(row_sum = sum(c_across(-id), na.rm = TRUE)) %>%
        mutate(across(-c(id, row_sum), ~./row_sum)) %>%
        select(-row_sum) %>% ungroup()
        
# Save results:
save(rna.seq.exp.profils, rpkm.rna.seq.counts, file = file.path(output_data_dir,"gene_expression.RData"))

write.table(rna.seq.exp.profils, file.path(output_data_dir, "RNA_Seq_RPKM_and_profiles.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

message("DONE")

