require(GeneFamilies)
options(mc.cores = getMcCores())

library(dotenv)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("USAGE: Rscript path/2/GeneFamilies/exec/loadRpkmRnaSeqCounts.R RPKM_counts_table.tsv")
message("Note, that the table in argument file 'RPKM_counts_table.tsv' is expected to be TAB-Delimited and have a header line:\n", 
    "id tissue rank expression variance")

input.args <- commandArgs(trailingOnly = TRUE)

#' Read RPKM normalized counts:
rpkm.rna.seq.counts <- read.table(input.args[[1]], sep = "\t", header = TRUE, stringsAsFactors = FALSE, 
    comment.char = "", quote = "", na.strings = "", colClasses = c(rep("character", 
        3), rep("numeric", "2")))
#' Transform into expression profiles:
genes <- sort(unique(rpkm.rna.seq.counts$id))
tissues <- sort(unique(rpkm.rna.seq.counts$tissue))

# All posible issued
all_tissues <- unique(rpkm.rna.seq.counts$tissue)

# Process every gene
process_gene <- function(x) {
    y <- rpkm.rna.seq.counts[which(rpkm.rna.seq.counts$id == x), ]
    
    # dataframe with all possible tissues
    x.df <- as.data.frame(matrix(NA, nrow = 1, ncol = length(all_tissues)))
    colnames(x.df) <- all_tissues
    rownames(x.df) <- NULL
    
    # Fill the tissues that exists for this gene
    tissue_names <- y$tissue
    expressions <- y$expression / sum(y$expression, na.rm = TRUE)
    
    for (i in seq_along(tissue_names)) {
        x.df[tissue_names[i]] <- expressions[i]
    }
    
    # Add gene column
    x.df$gene <- x
    x.df
}

# Combine results
rna.seq.exp.profils <- tryCatch({
    do.call("rbind", mclapply(genes, process_gene, mc.cores = parallel::detectCores()))  # Ajusta el número de núcleos según tus necesidades
}, error = function(e) {
    cat("Error al combinar los data.frames con rbind:\n")
    cat("Mensaje de error:", e$message, "\n")
    NULL  # Retorna NULL en caso de error
})

# #' Add matching gene-names _with_ expression variants:
# rna.seq.exp.profils$gene.exp.var <- as.character(unlist(mclapply(rna.seq.exp.profils$gene, 
#     function(x) {
#         x.exp.var <- names(all.cds)[grepl(paste("^", x, sep = ""), names(all.cds))]
#         if (length(x.exp.var) == 1) {
#             x.exp.var
#         } else {
#             warning("Found != 1 matching expression variants for '", x, "': ", paste(x.exp.var, 
#                 collapse = ", "), " !")
#             NA
#         }
#     })))

#' Save results:
save(rna.seq.exp.profils, rpkm.rna.seq.counts, file = file.path(output_data_dir, "RNA_Seq_RPKM_and_profiles.RData"))
write.table(rna.seq.exp.profils, file.path(output_data_dir, "RNA_Seq_RPKM_and_profiles.tsv"), 
    sep = "\t", row.names = FALSE, quote = FALSE)

message("DONE")
