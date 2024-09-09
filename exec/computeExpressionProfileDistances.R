require(GeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/GeneFamilies/exec/computeExpressionProfileDistances.R path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)
load("data/families.RData")
load("data/RNA_Seq_RPKM_and_profiles.RData")
load("data/orthologsTandems.RData")

expressionProfilesDists_Viv <- function(gene.accessions, expression.profiles = rna.seq.exp.profils, 
    expr.prof.gene.col = "gene", 
    tissues = setdiff(colnames(expression.profiles), c(expr.prof.gene.col)), 
    dist.method = "euclidean", per.tissue = FALSE) {
    
    # Encuentra los índices de los genes en gene.accessions
    inds <- which(expression.profiles[, expr.prof.gene.col] %in% gene.accessions)
    # print(gene.accessions)
    # print(inds)
    # Si se encuentran más de un gen, se calculan las distancias
    if (length(inds) > 1) {
        exp.profs <- expression.profiles[inds, ]
        rownames(exp.profs) <- exp.profs[, expr.prof.gene.col]
        
        # Si per.tissue es TRUE, calcula la distancia por tejido
        if (per.tissue) {
            setNames(lapply(tissues, function(tissue) dist(setNames(exp.profs[, 
                tissue], rownames(exp.profs)), method = dist.method)), 
                tissues)
        } else {
            # Calcula la distancia global
            dist(exp.profs[, tissues], method = dist.method)
        }
    } else {
        NA
    }
}


#' Compute for each group pairwise gene expression profile distances:
#' - Tandems:
# tandems.exp.prof.dists <- mclapply(tandems.lst, expressionProfilesDists)
# tandems.exp.prof.dists.tissue <- mclapply(tandems.lst, expressionProfilesDists, per.tissue = TRUE)
#' - Orthologs:
# orthologs.exp.prof.dists <- mclapply(orthologs.lst, expressionProfilesDists_Viv)
# orthologs.exp.prof.dists.tissue <- mclapply(orthologs.lst, expressionProfilesDists_Viv, 
#     per.tissue = TRUE)
# #' - Paralogs:
# print(paralogs.lst)
paralogs.exp.prof.dists <- mclapply(paralogs.lst, expressionProfilesDists_Viv)
paralogs.exp.prof.dists.tissue <- mclapply(paralogs.lst, expressionProfilesDists_Viv, 
    per.tissue = TRUE)
#' - Gene Families:
# print(head(families.df))
non.singleton.fams <- families.df$id[which(families.df$size > 1)]
# print(non.singleton.fams)
# print(families.lst)
# print(families.lst[non.singleton.fams])
families.exp.prof.dists <- mclapply(families.lst[non.singleton.fams], expressionProfilesDists_Viv)

families.exp.prof.dists.tissue <- mclapply(families.lst[non.singleton.fams], expressionProfilesDists_Viv, 
    per.tissue = TRUE)


#' Save results:
# save(tandems.exp.prof.dists, paralogs.exp.prof.dists.tissue, orthologs.exp.prof.dists, 
#     orthologs.exp.prof.dists.tissue, families.exp.prof.dists, families.exp.prof.dists.tissue, 
#     file = file.path(input.args[[1]], "data", "ExpressionProfileDistances.RData"))
save(paralogs.exp.prof.dists, paralogs.exp.prof.dists.tissue,families.exp.prof.dists, families.exp.prof.dists.tissue, 
    file = file.path(input.args[[1]], "data", "ExpressionProfileDistances.RData"))


message("DONE")
