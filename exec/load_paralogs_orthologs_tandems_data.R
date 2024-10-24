require(GeneFamilies)
library(dotenv)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("USAGE: Rscript exec/load_paralogs_orthologs_tandems_data.R <all_vs_all_file> <orthologs_file> <paralogs_file>")
message("input.args[[1]]: <all_vs_all_file> is a tabular Blast Output txt file, 
        resulted from running BLAST Tool on all coding sequences in an 'all vs all' approach.")
message("input.args[[2]]: <orthologs_file> is a TAB separated table or txt file")
message("input.args[[3]]: <paralogs_file> is a TAb separated table or txt file")

input.args <- commandArgs(trailingOnly = TRUE)

# Load data:
all.vs.all.sim <- fread(input.args[[1]], data.table = FALSE, header = FALSE, stringsAsFactors = FALSE, 
    sep = "\t", na.strings = "", colClasses = c(rep("character", 2), rep("numeric", 10)))

orthologs <- read.table(input.args[[1]], header = TRUE, sep = "\t", comment.char = "", 
    quote = "", na.strings = "", colClasses = rep("character", 1))
orths.nms <- paste("ortholog_cluster_", 1:nrow(orthologs), sep = "")
orthologs.lst <- setNames(mclapply(1:nrow(orthologs), function(x) {
    unlist(strsplit(orthologs[x, 1], ",\\s*"))  
}), orths.nms)
orthologs.genes <- unlist(orthologs.lst)

paralogs <- read.table(input.args[[2]], header = TRUE, sep = "\t", comment.char = "", 
    quote = "", na.strings = "", colClasses = rep("character", 2))
paralogs.nms <- unique(paralogs$Family)
paralogs.lst <- setNames(mclapply(paralogs.nms, function(x) paralogs[which(paralogs$Family == 
    x), "Gene"]), paralogs.nms)
paralogs.genes <- unlist(paralogs.lst)
paralogs.genes.nev <- removeExpressionVariant(paralogs.genes)

# tandems <- read.table(input.args[[3]], header = TRUE, sep = "\t", comment.char = "", 
#     quote = "", na.strings = "", colClasses = rep("character", 2))

# Save data:
save(orthologs, orthologs.lst,orthologs.genes,paralogs, paralogs.lst, paralogs.genes, paralogs.genes.nev, 
    file = file.path(output_data_dir, "orthologsTandems.RData"))

save(all.vs.all.sim, file = file.path(output_data_dir, "pairwiseSequenceSimilarities.RData"))

message("DONE")

# All vs All is missing ??