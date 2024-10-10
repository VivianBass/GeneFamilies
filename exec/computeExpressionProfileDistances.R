require(GeneFamilies)
options(mc.cores = getMcCores())
library(dotenv)
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("USAGE: Rscript path/2/GeneFamilies/exec/computeExpressionProfileDistances.R")

input.args <- commandArgs(trailingOnly = TRUE)
load(file.path(output_data_dir, "families.RData"))      # contiene el objeto de familias con columnas Family y Gene
load(file.path(output_data_dir, "orthologsTandems.RData")) # contiene paralogs.lst y orthologs.lst
load(file.path(output_data_dir, "RNA_Seq_RPKM_and_profiles.RData"))      # contiene las expresiones genéticas por tissue, con los genes en la última columna

test_expressionProfilesDists <- function(gene.accessions, expression.profiles = rna.seq.exp.profils,
 expr.prof.gene.col = "gene",
    tissues = setdiff(colnames(expression.profiles), c(expr.prof.gene.col)),
    dist.method = "euclidean", per.tissue = FALSE) {

    all_genes <- unlist(gene.accessions)
    # print(all_genes)
    inds <- which(expression.profiles$gene %in% all_genes)

    # print(length(inds))
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



orthologs.exp.prof.dists <- mclapply(orthologs.lst, test_expressionProfilesDists)
orthologs.exp.prof.dists.tissue <- mclapply(orthologs.lst, test_expressionProfilesDists,
     per.tissue = TRUE)
# #' - Paralogs:
# print(paralogs.lst)
paralogs.exp.prof.dists <- mclapply(paralogs.lst, test_expressionProfilesDists)
paralogs.exp.prof.dists.tissue <- mclapply(paralogs.lst, test_expressionProfilesDists,
   per.tissue = TRUE)
#' - Gene Families:
# print(head(families.df))
non.singleton.fams <- families.df$id[which(families.df$size > 1)]
# print(non.singleton.fams)
# print(families.lst)
# print(families.lst[non.singleton.fams])
families.exp.prof.dists <- mclapply(families.lst[non.singleton.fams], test_expressionProfilesDists)

families.exp.prof.dists.tissue <- mclapply(families.lst[non.singleton.fams], test_expressionProfilesDists,
    per.tissue = TRUE)


#' Save results:
# save(tandems.exp.prof.dists, paralogs.exp.prof.dists.tissue, orthologs.exp.prof.dists, 
#     orthologs.exp.prof.dists.tissue, families.exp.prof.dists, families.exp.prof.dists.tissue, 
#     file = file.path(input.args[[1]], "data", "ExpressionProfileDistances.RData"))
save(orthologs.exp.prof.dists,orthologs.exp.prof.dists.tissue,
paralogs.exp.prof.dists, paralogs.exp.prof.dists.tissue,families.exp.prof.dists, families.exp.prof.dists.tissue,
    file = file.path(output_data_dir, "ExpressionProfileDistances.RData"))


message("DONE")

