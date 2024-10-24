require(GeneFamilies)
options(mc.cores = getMcCores())

library(dotenv)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("USAGE:  Rscript exec/load_expression_data.R <RPKM_counts_table.tsv>")
message("Note:   <RPKM_counts_table.tsv> is expected to be TAB-Delimited and have a header line:\n", "id | tissue | rank")

input.args <- commandArgs(trailingOnly = TRUE)

# read RPKM normalized counts:
rpkm.rna.seq.counts <- read.table(input.args[[1]], sep = "\t", header = TRUE, stringsAsFactors = FALSE, 
    comment.char = "", quote = "", na.strings = "", colClasses = c(rep("character", 3), rep("numeric", "2")))

# Transform into expression profiles:
genes <- sort(unique(rpkm.rna.seq.counts$id))
tissues <- sort(unique(rpkm.rna.seq.counts$tissue))
all_tissues <- unique(rpkm.rna.seq.counts$tissue)

# Process every gene
process_gene <- function(x) {
    y <- rpkm.rna.seq.counts[which(rpkm.rna.seq.counts$id == x), ]
    
    # dataframe with all possible tissues
    x.df <- as.data.frame(matrix(NA, nrow = 1, ncol = length(all_tissues)))
    colnames(x.df) <- all_tissues
    rownames(x.df) <- NULL
    
    tissue_names <- y$tissue
    expressions <- y$expression / sum(y$expression, na.rm = TRUE)
    
    for (i in seq_along(tissue_names)) {
        x.df[tissue_names[i]] <- expressions[i]
    }
    
    x.df$gene <- x
    x.df
}

rna.seq.exp.profils <- tryCatch({
    do.call("rbind", mclapply(genes, process_gene, mc.cores = parallel::detectCores()))  
}, error = function(e) {
    cat("Error al combinar los data.frames con rbind:\n")
    cat("Mensaje de error:", e$message, "\n")
    NULL  
})

# Save results:
save(rna.seq.exp.profils, rpkm.rna.seq.counts, file = file.path(output_data_dir, "RNA_Seq_RPKM_and_profiles.RData"))
write.table(rna.seq.exp.profils, file.path(output_data_dir, "RNA_Seq_RPKM_and_profiles.tsv"), 
    sep = "\t", row.names = FALSE, quote = FALSE)

message("DONE")
