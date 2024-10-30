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
message("<RPKM_counts_table.tsv> Header : \n", "id | tissue | expression")

input.args <- commandArgs(trailingOnly = TRUE)

input.args[[1]] <- "experiments/test/RPKM.tsv"

# read RPKM counts:
rpkm.rna.seq.counts <- read.table(input.args[[1]], 
        sep = "\t", header = TRUE, fill = TRUE,
        stringsAsFactors = FALSE, comment.char = "", quote = "", na.strings = "") %>%
        select(id, tissue, expression) %>%
        mutate(expression = as.numeric(expression))

genes <- sort(unique(rpkm.rna.seq.counts$id))
tissues <- sort(unique(rpkm.rna.seq.counts$tissue))

# compute expression profiles for each gene and normalize them
rna.seq.exp.profils <- do.call("rbind", mclapply(genes, function(x) {
    y <- rpkm.rna.seq.counts[which(rpkm.rna.seq.counts$id == x), ]
    x.df <- as.data.frame(t(setNames(y[, "expression"]/sum(y[, "expression"], 
    na.rm = TRUE),  y$tissue)), stringsAsFactors = FALSE)
    x.df$gene <- x
    x.df
}))

# Save results:
save(rna.seq.exp.profils, rpkm.rna.seq.counts, file = file.path(output_data_dir,"gene_expression.RData"))

write.table(rna.seq.exp.profils, file.path(output_data_dir, "RNA_Seq_RPKM_and_profiles.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

message("DONE")
