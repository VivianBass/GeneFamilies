require(GeneFamilies)
options(mc.cores = getMcCores())
library(dotenv)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("USAGE: Rscript exec/compute_exp.prof.dists.R")

input.args <- commandArgs(trailingOnly = TRUE)
load(file.path(output_data_dir, "gene_families.RData"))                       
load(file.path(output_data_dir, "gene_groups.RData"))              
load(file.path(output_data_dir, "gene_expression.RData"))     

# Function for ...
expressionProfilesDists <- function(gene.accessions, expression.profiles = rna.seq.exp.profils,
    expr.prof.gene.col = "gene",
    tissues = setdiff(colnames(expression.profiles), c(expr.prof.gene.col)),
    dist.method = "euclidean", per.tissue = FALSE) {

    all_genes <- unlist(gene.accessions)
    inds <- which(expression.profiles$gene %in% all_genes)

    if (length(inds) > 1) {

        exp.profs <- expression.profiles[inds, ]
        exp.profs <- as.data.frame(exp.profs)
        rownames(exp.profs) <- exp.profs$gene

        if (per.tissue) {
            setNames(mclapply(tissues, function(tissue) dist(setNames(exp.profs[,
                tissue], rownames(exp.profs)), method = dist.method)),
                tissues)
        } else dist(exp.profs[, tissues], method = dist.method)

    } else NA
}

# Orthologs:
orthologs.exp.prof.dists <- mclapply(orthologs.lst, expressionProfilesDists)
orthologs.exp.prof.dists.tissue <- mclapply(orthologs.lst, expressionProfilesDists,
                                            per.tissue = TRUE)

# Paralogs:
paralogs.exp.prof.dists <- mclapply(paralogs.lst, expressionProfilesDists)
paralogs.exp.prof.dists.tissue <- mclapply(paralogs.lst, expressionProfilesDists,
                                            per.tissue = TRUE)

# Gene-Families:
non.singleton.fams <- families.df$id[which(families.df$size > 1)]
families.exp.prof.dists <- mclapply(families.lst[non.singleton.fams], expressionProfilesDists)
families.exp.prof.dists.tissue <- mclapply(families.lst[non.singleton.fams], expressionProfilesDists,
                                            per.tissue = TRUE)


save(orthologs.exp.prof.dists,orthologs.exp.prof.dists.tissue, 
     paralogs.exp.prof.dists, paralogs.exp.prof.dists.tissue,
     families.exp.prof.dists, families.exp.prof.dists.tissue,
     file = file.path(output_data_dir, "exp.prof.dists.RData"))

message("DONE")

