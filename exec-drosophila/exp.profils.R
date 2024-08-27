

rpkm.rna.seq.counts <- read.table("p.df.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                                  comment.char = "", quote = "", na.strings = "", 
                                  colClasses = c(rep("character", 3), rep("numeric", 2)))

# Get unique genes and tissues
genes <- sort(unique(rpkm.rna.seq.counts$id))
tissues <- sort(unique(rpkm.rna.seq.counts$tissue))

# Create a function to process each gene
process_gene <- function(x) {
  y <- rpkm.rna.seq.counts[which(rpkm.rna.seq.counts$id == x), ]
  
  # Create a data frame with all tissues, filling missing values with 0
  x.df <- data.frame(matrix(0, nrow = 1, ncol = length(tissues)))
  colnames(x.df) <- tissues
  
  # Fill in the available data
  x.df[1, y$tissue] <- y$expression / sum(y$expression, na.rm = TRUE)
  
  x.df$gene <- x
  return(x.df)
}

# Apply the function to all genes
library(parallel)
rna.seq.exp.profils <- do.call("rbind", mclapply(genes, process_gene))


genes_no_exp <- rpkm.rna.seq.counts[rpkm.rna.seq.counts$expression == 0 | is.na(rpkm.rna.seq.counts$expression), "id"]
print(unique(genes_no_exp))