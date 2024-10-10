require(GeneFamilies)
options(mc.cores = getMcCores())
library(dotenv)
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")


message("USAGE: Rscript path/2/GeneFamilies/exec/investigateDistributionsOfExpressionProfileDistances.R path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)

load(file.path(output_data_dir, "orthologsTandems.RData")) # contiene paralogs.lst y orthologs.lst
load(file.path(output_data_dir, "families.RData"))      # contiene el objeto de familias con columnas Family y Gene
load(file.path(output_data_dir,"pairwiseSequenceSimilarities.RData"))
load(file.path(output_data_dir,"ExpressionProfileDistances.RData"))


classSpecificExpressionProfileDists_test <- function(e.p.d, gene.classes, fam.name = NA, 
    stats = c("max", "median", "mean", "maxMinusMin")) {
    epd.m <- as.matrix(e.p.d)
    genes <- colnames(epd.m)
    genes.classified <- lapply(gene.classes, function(x) intersect(genes, 
        x))
    # print(genes.classified)
    res <- lapply(genes.classified, function(x) as.dist(epd.m[rownames(epd.m) %in% 
        x, colnames(epd.m) %in% x]))
    res$stats <- Reduce(rbind, lapply(stats, function(i.stat) {
        Reduce(rbind, lapply(names(gene.classes), function(k.class) {
            genes.dist <- res[[k.class]]
            dist.stat <- if (length(genes.dist) > 0 && !all(is.na(genes.dist) | 
                is.nan(genes.dist) | is.infinite(genes.dist))) {
                eval(call(i.stat, unlist(res[[k.class]]), na.rm = TRUE))
            } else NA
            data.frame(Family = fam.name, Statistic = i.stat, Gene.Class = k.class, 
                Value = dist.stat, stringsAsFactors = FALSE)
        }))
    }))
    res
}



#' Compute statistics about Expression Profile Distances within different
#' classes of gene groups. In this, consider orthologs to be the 'background'.
gene.classes <- list(Orthologs = sub("\\.\\d+$", "", orthologs.genes), `Paralogs` = sub("\\.\\d+$", 
    "", sub("\\.\\d+$", "", paralogs.genes)))
# print(gene.classes)
families_with_paralogs <- vector()

# Itera sobre cada familia de genes
for (family_name in names(families.exp.prof.dists)) {
    # Obtén los genes de la familia actual
    family_genes <- colnames(families.exp.prof.dists[[family_name]])
    
    # Busca coincidencias con los parálogos
    paralogs_in_family <- intersect(family_genes, gene.classes$Paralogs)
    
    # Si hay parálogos en la familia, agrega el nombre de la familia a la lista
    if (length(paralogs_in_family) > 0) {
        families_with_paralogs <- c(families_with_paralogs, family_name)
        # Opcionalmente, imprimir los parálogos encontrados
        print(paste("Familia:", family_name, "contiene parálogos:", paste(paralogs_in_family, collapse = ", ")))
    }
}

# Imprimir las familias que contienen parálogos
print("Familias que contienen parálogos:")
print(families_with_paralogs)

#' - Tandems:
# tandems.exp.prof.dists.orth.dist <- setNames(mclapply(names(tandems.exp.prof.dists), 
#     function(x) classSpecificExpressionProfileDists(tandems.exp.prof.dists[[x]], 
#         gene.classes = gene.classes, fam.name = x)), names(tandems.exp.prof.dists))
# tandems.exp.prof.dists.orth.dist.df <- Reduce(rbind, mclapply(names(tandems.exp.prof.dists.orth.dist), 
#     function(x) tandems.exp.prof.dists.orth.dist[[x]][["stats"]]))
#' - Families:
# families.exp.prof.dists.orth.dist <- setNames(mclapply(names(families.exp.prof.dists), 
#     function(x) classSpecificExpressionProfileDists_test(families.exp.prof.dists[[x]], 
#         gene.classes = gene.classes, fam.name = x)), names(families.exp.prof.dists))
# families.exp.prof.dists.orth.dist.df <- Reduce(rbind, mclapply(names(families.exp.prof.dists.orth.dist), 
#     function(x) families.exp.prof.dists.orth.dist[[x]][["stats"]]))
#' - Orthologs:
#' (Of course no distinction can be made between Non-Orthologs and Orthologs
#' within _Orthologs_ themselves, but none the less we need the statistics for
#' each cluster of orthologs.)
# orthologs.exp.prof.dists.stats <- setNames(mclapply(names(orthologs.exp.prof.dists), 
#     function(x) classSpecificExpressionProfileDists_test(orthologs.exp.prof.dists[[x]], 
#         gene.classes = gene.classes["Orthologs"], fam.name = x)), names(orthologs.exp.prof.dists))
# orthologs.exp.prof.dists.stats.df <- Reduce(rbind, mclapply(names(orthologs.exp.prof.dists.stats), 
#     function(x) orthologs.exp.prof.dists.stats[[x]][["stats"]]))


# paralogs.exp.prof.dists.stats <- setNames(mclapply(names(paralogs.exp.prof.dists), 
#     function(x) classSpecificExpressionProfileDists_test(paralogs.exp.prof.dists[[x]], 
#         gene.classes = gene.classes["Paralogs"], fam.name = x)), names(paralogs.exp.prof.dists))
# paralogs.exp.prof.dists.stats.df <- Reduce(rbind, mclapply(names(paralogs.exp.prof.dists.stats), 
#     function(x) paralogs.exp.prof.dists.stats[[x]][["stats"]]))


#' Save results:
# save(paralogs.exp.prof.dists.stats.df, families.exp.prof.dists.orth.dist, 
#     families.exp.prof.dists.orth.dist.df, orthologs.exp.prof.dists.stats.df, file = file.path(output_data_dir, "ExpressionProfileDistanceDistributions.RData"))


message("DONE")
