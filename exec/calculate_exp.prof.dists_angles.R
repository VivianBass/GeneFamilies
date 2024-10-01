
# calculate angles between expression profiles and diagonal

# Functions used in this script from expression_functions.R:

# cosAngleVec (from calculate_expression_profile_dists_function.R)
# cosDiag (from calculate_expression_profile_dists_function.R) 
# euclNorm (from calculate_expression_profile_dists_function.R)
# hasInvalidVecSpaceComponents (from calculate_expression_profile_dists_function.R)
# addalpha (from expression_functions.R)

# plots generated:
# "expressionAngleToDiagonalBoxplot.pdf"
# "relativeExpressionVersatilityBoxplot.pdf"

# further Functions from expression_functions.R for further Analysis:

# scalarProjection (from calculate_expression_profile_dists_function.R)
# vectorProjection (from calculate_expression_profile_dists_function.R)
# statVectorCloud (from calculate_expression_profile_dists_function.R)
# distVectorClouds (from calculate_expression_profile_dists_function.R)
# naAsZero (from calculate_expression_profile_dists_function.R)
# exprVecSpaceEvolutionAfterDupl

# -------------------------------------------------------------------------------

# Inputs used in this script:
rpkm.expr.profiles.df <- df14_P_dmel
paralogs.lst <- vv_df3_para.lst
orthologs.lst <- vv_df3_ortho.lst

# -------------------------------------------------------------------------------

library(parallel)
source("R/expression_funks.R")

#' Define expression vector space:
expr.cols <- c("Whole_Body", "Muscle", "Gut", "Fat_Body")
n.dims <- length(expr.cols)

# calculate exp.prof.dists angles to diagonal for paralogs 

paralog.genes <- unlist(paralogs.lst)
para.expr <- intersect(paralog.genes, rpkm.expr.profiles.df$Gene)

paralog.expr.angle.diag.df <- data.frame(gene = para.expr, angle.diag = as.numeric(mclapply(para.expr, 
    function(x) {
        cosDiag(rpkm.expr.profiles.df[which(rpkm.expr.profiles.df$Gene == 
            x), expr.cols])/sqrt(2)
    })), stringsAsFactors = FALSE)


# calculate exp.prof.dists angles to diagonal for orthologs 

orthologs.genes <- unlist(orthologs.lst)
orths.expr <- intersect(orthologs.genes, rpkm.expr.profiles.df$Gene)

orths.expr.angle.diag.df <- data.frame(Gene = orths.expr, angle.diag = as.numeric(mclapply(orths.expr, 
    function(x) {
        cosDiag(rpkm.expr.profiles.df[which(rpkm.expr.profiles.df$Gene == 
            x), expr.cols])/sqrt(2)
    })), stringsAsFactors = FALSE)

# -------------------------------------------------------------------------------

# merge paralog.expr.angle.diag.df and orths.expr.angle.diag.df and plot results:

p.lst <- list(paralog = paralog.expr.angle.diag.df, ortholog = orths.expr.angle.diag.df)
p.df <- Reduce(rbind, mclapply(names(p.lst), function(gene.type) {
    data.frame(gene.type = gene.type, angle.diag = p.lst[[gene.type]]$angle.diag, 
        stringsAsFactors = FALSE)
}))

library(RColorBrewer)

plot.df <- p.df[!is.nan(p.df$angle.diag), ]
plot.df$gene.type <- factor(plot.df$gene.type, levels = c("ortholog", "paralog"))
pdf(file.path("expressionAngleToDiagonalBoxplot.pdf"))
colors <- brewer.pal(length(p.lst), "Dark2")
pushed.colors <- append(colors[3], colors[1:2])
boxplot(angle.diag ~ gene.type, data = plot.df, xlab = "Type of Gene", 
    ylab = "relative tissue specificity", border = pushed.colors, col = addAlpha(pushed.colors))
dev.off()


plot.df$rel.vers <- 1 - plot.df$angle.diag
pdf(file.path("relativeExpressionVersatilityBoxplot.pdf"))
boxplot(rel.vers ~ gene.type, data = plot.df, xlab = "Type of Gene", ylab = "relative tissue versatility", 
    border = pushed.colors, col = addAlpha(pushed.colors), outline = FALSE)
dev.off()

# plots generated:

# "expressionAngleToDiagonalBoxplot.pdf"
# "relativeExpressionVersatilityBoxplot.pdf"

# -------------------------------------------------------------------------------



tandems.list
orthologs.genes

#' Define Duplicated Gene Sets with Orthologs:
t.i <- which(as.logical(lapply(tandems.lst, function(t.c) any(removeExpressionVariant(t.c) %in% 
    orthologs.genes))))
tands.w.orths <- lapply(tandems.lst[t.i], function(t.c) removeExpressionVariant(t.c))
d.i <- which(as.logical(lapply(families.lst[families.exp.cnt.unq], function(f.e) any(removeExpressionVariant(f.e) %in% 
    orthologs.genes))))
dupl.w.orths <- lapply(families.lst[families.exp.cnt.unq][d.i], function(f.e) removeExpressionVariant(f.e))


tands.w.orths
rpkm.expr.profiles.df

#' Investigate what happened in each gene group after duplication:
#' - Comparing tandems and orthologs
tands.w.orths.angles.df <- Reduce(rbind, mclapply(names(tands.w.orths), 
    function(fam.name) {
        genes <- intersect(tands.w.orths[[fam.name]], rpkm.expr.profiles.df$gene)
        exprVecSpaceEvolutionAfterDupl(genes, fam.name, tand.classifier)
    }))

tands.w.orths.angles.df
dupl.w.orths.angles.df

#' Plot the above:
p.lst <- list(tandem = tands.w.orths.angles.df[which(tands.w.orths.angles.df$tandem.dist.vec.clouds > 
    0), ], duplicated = dupl.w.orths.angles.df[which(dupl.w.orths.angles.df$duplicated.dist.vec.clouds > 
    0), ])

p.df <- Reduce(rbind, lapply(names(p.lst), function(x) {
    tryCatch({
        y <- p.lst[[x]]
        data.frame(group.type = x, diff.tissue.versatility = (y[, paste(x, 
            ".tiss.vers", sep = "")] - y$mean.base.tiss.vers), angle.between.orth.2.diag.vecs = y[, 
            paste(x, ".tiss.change", sep = "")], stringsAsFactors = FALSE)
    }, error = function(e) browser())
}))
"meanTissueVersatilityDiffsAfterDuplicationBoxplot.pdf"
"meanTissueVersatilityDiffsAfterDuplicationHistograms.pdf"

pdf(file.path(input.args[[1]], "inst", "meanTissueVersatilityDiffsAfterDuplicationBoxplot.pdf"))
boxplot(diff.tissue.versatility ~ group.type, data = p.df, col = addAlpha(colors), 
    outline = FALSE, border = colors, xlab = "Type of gene group", ylab = "Mean difference of tissue versatility after duplication", 
    pch = "-")
stripchart(diff.tissue.versatility ~ group.type, vertical = TRUE, data = p.df, 
    method = "jitter", add = TRUE, pch = ".", col = colors)
dev.off()


pdf(file.path(input.args[[1]], "inst", "meanTissueVersatilityDiffsAfterDuplicationHistograms.pdf"))
old.par <- par(mfrow = c(2, 1))
xlim <- c(min(p.df$diff.tissue.versatility, na.rm = TRUE), max(p.df$diff.tissue.versatility, 
    na.rm = TRUE))
hist(p.df[which(p.df$group.type == "duplicated"), "diff.tissue.versatility"], 
    breaks = 30, main = "duplicated", col = addAlpha(colors[[1]]), border = colors[[1]], 
    xlab = "", xlim = xlim)
hist(p.df[which(p.df$group.type == "tandem"), "diff.tissue.versatility"], 
    breaks = 30, main = "tandem", col = addAlpha(colors[[2]]), border = colors[[2]], 
    xlab = "Mean difference of tissue versatility after duplication", xlim = xlim)
dev.off()
par(old.par)

"meanTissueVersatilityDiffsAfterDuplicationBoxplot.pdf"
"meanTissueVersatilityDiffsAfterDuplicationHistograms.pdf"