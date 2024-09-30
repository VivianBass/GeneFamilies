



load("x/data/orthologs/orthologs_standart.RData")


source("x/Section-2/family_funks.R")



orths.nms <- paste("ortholog_cluster_", 1:nrow(ortho), sep = "")
ortho.lst <- setNames(mclapply(1:nrow(ortho), function(x) unlist(ortho[x, ])), orths.nms)
ortho.genes <- unlist(ortho.lst)




#' - Orthologs:
orths.nms <- paste("ortholog_cluster_", 1:nrow(orthologs), sep = "")
orthologs.lst <- setNames(mclapply(1:nrow(orthologs), function(x) unlist(orthologs[x, ])), orths.nms)
orthologs.genes <- unlist(orthologs.lst)


# Wie am besten Paralogs bestimmen ?




#' - Paralogs V-1  XXX :
paralogs.nms <- unique(paralogs$Family)
paralogs.lst <- setNames(mclapply(paralogs.nms, function(x) paralogs[which(paralogs$Family == 
    x), "Gene"]), paralogs.nms)
paralogs.genes <- unlist(paralogs.lst)
paralogs.genes.nev <- removeExpressionVariant(paralogs.genes)




#' - Paralogs V-1:
paras.nms <- paste("paralog_cluster_", 1:nrow(paralogs), sep = "")
paralogs.lst <- setNames(mclapply(1:nrow(paralogs), function(x) unlist(paralogs[x, ])), paras.nms)
paralogs.genes <- unlist(paralogs.lst)


save(orthologs, orthologs.lst , file = file.path("orthologs.RData"))
save(paralogs, paralogs.lst , file = file.path("paralogs.RData"))