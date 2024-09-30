require(GeneFamilies)
options(mc.cores = getMcCores())

library(parallel)

message("USAGE: Rscript path/2/GeneFamilies/exec/computeExpressionProfileDistances.R path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)

# rna.seq.exp.profils
# paralogs.lst, orthologs.lst, families.lst


load("new_data_X.RData")
load("ortholog_data_1.RData")
load("data/paralogs/paralogs.RData")
load("paralogs_filtered.RData")

source("data/R/expression_funks.R")



expressionProfilesDists <- function(gene.accessions, expression.profiles = rna.seq.exp.profils, 
    expr.prof.gene.col = "Gene", expr.prof.gene.variant.col = "gene.exp.var", 
    tissues = setdiff(colnames(expression.profiles), c(expr.prof.gene.col, 
        expr.prof.gene.variant.col)), dist.method = "euclidean", per.tissue = FALSE) {
    inds <- which(expression.profiles[, expr.prof.gene.col] %in% gene.accessions | 
        expression.profiles[, expr.prof.gene.variant.col] %in% gene.accessions)
    if (length(inds) > 1) {
        exp.profs <- expression.profiles[inds, ]
        rownames(exp.profs) <- exp.profs[, expr.prof.gene.col]
        if (per.tissue) {
            setNames(lapply(tissues, function(tissue) dist(setNames(exp.profs[, 
                tissue], rownames(exp.profs)), method = dist.method)), 
                tissues)
        } else dist(exp.profs[, tissues], method = dist.method)
    } else NA
}


orthologs.lst <- vv_df3_ortho.lst
rna.seq.exp.profils <- df14_P_dmel


expressionProfilesDists <- function(orthologs.lst, rna.seq.exp.profils, per.tissue = FALSE) {

    tissues = setdiff(colnames(rna.seq.exp.profils), "Gene")
    dist.method = "euclidean"

    all_genes <- unlist(orthologs.lst)
    inds <- which(rna.seq.exp.profils$Gene %in% all_genes)

    if (length(inds) > 1) {

        exp.profs <- rna.seq.exp.profils[inds, ]
        exp.profs <- as.data.frame(exp.profs)
        rownames(exp.profs) <- exp.profs$Gene

        if (per.tissue) {
            setNames(mclapply(tissues, function(tissue) dist(setNames(exp.profs[, 
                tissue], rownames(exp.profs)), method = dist.method)), 
                tissues)
        } else dist(exp.profs[, tissues], method = dist.method)

    } else NA
}

library(parallel)

orthologs.exp.prof.dists <- mclapply(orthologs.lst, function(x) expressionProfilesDists(x, rna.seq.exp.profils, per.tissue = FALSE))

orthologs.exp.prof.dists.tissue <- mclapply(orthologs.lst, function(x) expressionProfilesDists(x, rna.seq.exp.profils, per.tissue = TRUE))

save(orthologs.exp.prof.dists, orthologs.exp.prof.dists.tissue,  
    file = file.path("OrthologsExpressionProfileDistances_P_dmel.RData"))



paralogs.lst <- vv_df3_para.lst
rna.seq.exp.profils <- df14_P_dmel

expressionProfilesDists <- function(paralogs.lst, rna.seq.exp.profils, per.tissue = FALSE) {

    tissues = setdiff(colnames(rna.seq.exp.profils), "Gene")
    dist.method = "euclidean"

    all_genes <- unlist(paralogs.lst)
    inds <- which(rna.seq.exp.profils$Gene %in% all_genes)

    if (length(inds) > 1) {

        exp.profs <- rna.seq.exp.profils[inds, ]
        exp.profs <- as.data.frame(exp.profs)
        rownames(exp.profs) <- exp.profs$Gene

        if (per.tissue) {
            setNames(mclapply(tissues, function(tissue) dist(setNames(exp.profs[, 
                tissue], rownames(exp.profs)), method = dist.method)), 
                tissues)
        } else dist(exp.profs[, tissues], method = dist.method)

    } else NA
}

library(parallel)

paralogs.exp.prof.dists <- mclapply(paralogs.lst, function(x) expressionProfilesDists(x, rna.seq.exp.profils, per.tissue = FALSE))

paralogs.exp.prof.dists.tissue <- mclapply(paralogs.lst, function(x) expressionProfilesDists(x, rna.seq.exp.profils, per.tissue = TRUE))

save(paralogs.exp.prof.dists, paralogs.exp.prof.dists.tissue,  
    file = file.path("ParalogsExpressionProfileDistances_P_dmel.RData"))




families.df


expressionProfilesDists <- function(families.lst, rna.seq.exp.profils, per.tissue = FALSE) {

    tissues = setdiff(colnames(rna.seq.exp.profils), "Gene")
    dist.method = "euclidean"

    all_genes <- unlist(families.lst)
    inds <- which(rna.seq.exp.profils$Gene %in% all_genes)

    if (length(inds) > 1) {

        exp.profs <- rna.seq.exp.profils[inds, ]
        exp.profs <- as.data.frame(exp.profs)
        rownames(exp.profs) <- exp.profs$Gene

        if (per.tissue) {
            setNames(mclapply(tissues, function(tissue) dist(setNames(exp.profs[, 
                tissue], rownames(exp.profs)), method = dist.method)), 
                tissues)
        } else dist(exp.profs[, tissues], method = dist.method)

    } else NA
}

non.singleton.fams <- families.df$id[which(families.df$size > 1)]

families.exp.prof.dists <- mclapply(families.lst[non.singleton.fams], function(x) expressionProfilesDists(x, rna.seq.exp.profils, per.tissue = TRUE))

families.exp.prof.dists.tissue <- mclapply(families.lst[non.singleton.fams], function(x) expressionProfilesDists(x, rna.seq.exp.profils, per.tissue = TRUE))


save(families.exp.prof.dists, families.exp.prof.dists.tissue,  
    file = file.path("FamiliesExpressionProfileDistances_R.RData"))




save(paralogs.exp.prof.dists, paralogs.exp.prof.dists.tissue, orthologs.exp.prof.dists, 
    orthologs.exp.prof.dists.tissue, families.exp.prof.dists, families.exp.prof.dists.tissue, 
    file = file.path(input.args[[1]], "data", "ExpressionProfileDistances.RData"))