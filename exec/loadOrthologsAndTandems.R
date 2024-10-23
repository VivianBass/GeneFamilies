require(GeneFamilies)
library(dotenv)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")


message("USAGE: Rscript path/2/GeneFamilies/exec/loadOrthologsAndParalogs.R eight_brassicaceae_orthologs.txt eight_brassicaceae_tandems.txt")
message("EXPEXTED FORMATS:")
message("- all_vs_all_tabular_blast_out.txt is the result of calling BLAT or BLAST on all coding sequences in an 'all vs all' approach. The result file is expected to be in the tabular Blast output format ")
message("- eight_brassicaceae_orthologs.txt is a TAB separated table holding orthologous gene clusters. The table is expected to have the following header:\naly\tath\tcru\tchi\taet\tbra\tesa\ttpa\n")
message("- eight_brassicaceae_tandems.txt is a TAb separated table with header 'Family Gene' and rows in the form of 'tandem_cluster_1 CARHR1234.1'.")

input.args <- commandArgs(trailingOnly = TRUE)

input.args[[1]] <- "experiments/RPKM_flybase/orthologs.csv"
input.args[[2]] <- "experiments/RPKM_flybase/paralogs.csv"

#' Load data:
#all.vs.all.sim <- fread(input.args[[1]], data.table = FALSE, header = FALSE, stringsAsFactors = FALSE, 
#    sep = "\t", na.strings = "", colClasses = c(rep("character", 2), rep("numeric", 
#        10)))
orthologs <- read.table(input.args[[1]], header = TRUE, sep = "\t", comment.char = "", 
    quote = "", na.strings = "", colClasses = rep("character", 1))

orths.nms <- paste("ortholog_cluster_", 1:nrow(orthologs), sep = "")
# Separar los genes por comas y generar una lista donde cada entrada es un cluster de genes
orthologs.lst <- setNames(mclapply(1:nrow(orthologs), function(x) {
    # Usamos strsplit para separar los genes por comas y espacios
    unlist(strsplit(orthologs[x, 1], ",\\s*"))  # Asumiendo que los genes están en la segunda columna
}), orths.nms)

# Crear un vector unificado con todos los genes
orthologs.genes <- unlist(orthologs.lst)

# tandems <- read.table(input.args[[3]], header = TRUE, sep = "\t", comment.char = "", 
#     quote = "", na.strings = "", colClasses = rep("character", 2))
paralogs <- read.table(input.args[[2]], header = TRUE, sep = "\t", comment.char = "", 
    quote = "", na.strings = "", colClasses = rep("character", 2))
paralogs.nms <- unique(paralogs$Family)
paralogs.lst <- setNames(mclapply(paralogs.nms, function(x) paralogs[which(paralogs$Family == 
    x), "Gene"]), paralogs.nms)
paralogs.genes <- unlist(paralogs.lst)
paralogs.genes.nev <- removeExpressionVariant(paralogs.genes)
#' Save loaded data:
save(orthologs, orthologs.lst,orthologs.genes,paralogs, paralogs.lst, paralogs.genes, paralogs.genes.nev, file = file.path(output_data_dir, "orthologsTandems.RData"))
#save(all.vs.all.sim, file = file.path(output_data_dir, "pairwiseSequenceSimilarities.RData"))

message("DONE")
