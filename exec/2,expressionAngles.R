require(GeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/GeneFamilies/exec/expressionAngles.R path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' **************************************
#' Analysis of expression vector space. *
#' **************************************


hasInvalidVecSpaceComponents <- function(x) {
    UseMethod("hasInvalidVecSpaceComponents", x)
}

hasInvalidVecSpaceComponents.numeric <- function(x) {
    any(is.na(x) | is.nan(x) | is.null(x) | is.infinite(x))
}

hasInvalidVecSpaceComponents.matrix <- function(x) {
    any(as.logical(apply(x, 1, hasInvalidVecSpaceComponents)))
}

hasInvalidVecSpaceComponents.data.frame <- function(x) {
    any(as.logical(apply(x, 1, hasInvalidVecSpaceComponents)))
}


----------------------------

cosAngleVec <- function(a, b) {
    if (length(a) != length(b)) 
        stop("'cosAngleVec(a, b)': Argument vectors 'a' and 'b' are not of identical dimensions.")
    if (hasInvalidVecSpaceComponents(rbind(a, b))) 
        return(NaN)
    sum(a * b)/(euclNorm(a) * euclNorm(b))
}

cosDiag <- function(x, d.v = rep(1, length(x))) {
    cosAngleVec(x, d.v)
}

euclNorm <- function(x) {
    if (hasInvalidVecSpaceComponents(x)) 
        return(NaN)
    sqrt(sum(x^2))
}
-------------------------------------------------------------------------------

# Infer angle to diagonal as measure of tissue specificity:
# --> cosDiag Function to generate Angles

#' Define expression vector space:
expr.cols <- c("Whole_Body", "Muscle", "Gut", "Fat_Body")
n.dims <- length(expr.cols)

# paste(expr.cols, collapse = ", "))

rpkm.expr.profiles.df <- df14_P_dmel
paralogs.lst <- vv_df3_para.lst
orthologs.lst <- vv_df3_ortho.lst

# Paralogs
paralog.genes <- unlist(paralogs.lst)
para.expr <- intersect(paralog.genes, rpkm.expr.profiles.df$Gene)

library(parallel)

paralog.expr.angle.diag.df <- data.frame(gene = para.expr, angle.diag = as.numeric(mclapply(para.expr, 
    function(x) {
        cosDiag(rpkm.expr.profiles.df[which(rpkm.expr.profiles.df$Gene == 
            x), expr.cols])/sqrt(2)
    })), stringsAsFactors = FALSE)



# Orthhologs
orthologs.genes <- unlist(orthologs.lst)
orths.expr <- intersect(orthologs.genes, rpkm.expr.profiles.df$Gene)

orths.expr.angle.diag.df <- data.frame(Gene = orths.expr, angle.diag = as.numeric(mclapply(orths.expr, 
    function(x) {
        cosDiag(rpkm.expr.profiles.df[which(rpkm.expr.profiles.df$Gene == 
            x), expr.cols])/sqrt(2)
    })), stringsAsFactors = FALSE)


# - merge to combined df
paralog.expr.angle.diag.df
orths.expr.angle.diag.df
df_para_ortho_angles

#' Plot results:
p.lst <- list(paralog = paralog.expr.angle.diag.df, ortholog = orths.expr.angle.diag.df)
p.df <- Reduce(rbind, mclapply(names(p.lst), function(gene.type) {
    data.frame(gene.type = gene.type, angle.diag = p.lst[[gene.type]]$angle.diag, 
        stringsAsFactors = FALSE)
}))

- plot
 "expressionAngleToDiagonalBoxplot.pdf"
 "relativeExpressionVersatilityBoxplot.pdf"


addAlpha <- function(col, alpha = 0.1) {
    apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], 
        alpha = alpha))
}

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


--------------------------------------------------------------------------

hasInvalidVecSpaceComponents.data.frame <- function(x) {
    any(as.logical(apply(x, 1, hasInvalidVecSpaceComponents)))
}

scalarProjection <- function(a, b) {
    if (length(a) != length(b)) 
        stop("'scalarProjection(a, b)': Argument vectors 'a' and 'b' are not of identical dimensions.")
    if (hasInvalidVecSpaceComponents(rbind(a, b))) 
        return(NaN)
    
    euclNorm(a) * cosAngleVec(a, b)
}

vectorProjection <- function(a, b) {
    if (length(a) != length(b)) 
        stop("'vectorProjection(a, b)': Argument vectors 'a' and 'b' are not of identical dimensions.")
    if (hasInvalidVecSpaceComponents(rbind(a, b))) 
        return(NaN)
    scalarProjection(a, b) * b/euclNorm(b)
}


statVectorCloud <- function(vecs.df, stat = getOption("GeneFamilies.vector.cloud.stat", 
    mean), deviance.funk = getOption("GeneFamilies.vector.cloud.deviance", 
    sd)) {
    if ((class(vecs.df) == "data.frame" || class(vecs.df) == "matrix") && 
        nrow(vecs.df) > 0) {
        i <- which(!as.logical(apply(vecs.df, 1, hasInvalidVecSpaceComponents)))
        if (length(i) == 0) 
            return(NA)
        v.df <- vecs.df[i, , drop = FALSE]
        stat.vec <- as.numeric(apply(v.df, 2, stat))
        dev.vec <- as.numeric(apply(v.df, 2, deviance.funk))
        orth.on.diag.2.stat.vec <- stat.vec - vectorProjection(stat.vec, 
            rep(1, length(stat.vec)))
        list(stat.vec = stat.vec, deviance.vec = dev.vec, orth.on.diag.2.stat.vec = orth.on.diag.2.stat.vec)
    }
}

naAsZero <- function(x) {
    if (is.na(x)) 
        0 else x
}

distVectorClouds <- function(cl.a, cl.b) {
    if (is.null(cl.a) || is.null(cl.b)) 
        return(NA)
    cl.diff <- cl.a$stat.vec - cl.b$stat.vec
    dev.proj.a <- naAsZero(scalarProjection(cl.a$deviance.vec, cl.diff))
    dev.proj.b <- naAsZero(scalarProjection(cl.b$deviance.vec, cl.diff))
    euclNorm(cl.diff) - (abs(dev.proj.a) + abs(dev.proj.b))
}

exprVecSpaceEvolutionAfterDupl <- function(genes, family.name, classifier.lst, 
    base.class = "ortholog", expr.vecs = rpkm.expr.profiles.df, gene.col = "Gene", 
    vec.space.axes = c("Whole_Body", "Muscle", "Gut", "Fat_Body")) {
    res <- NULL
    if (length(genes) > 0) {
        genes.classed <- lapply(classifier.lst, function(x) intersect(intersect(genes, 
            x), expr.vecs[, gene.col]))
        g.c.i <- sapply(genes.classed, function(x) length(x) > 0)
        evolved.classes <- setdiff(names(genes.classed), base.class)
        # If base class and at least one other class has genes with expression
        # values:
        if (length(genes.classed[[base.class]]) > 0 && any(g.c.i[evolved.classes])) {
            base.class.cloud <- statVectorCloud(expr.vecs[which(expr.vecs[, 
                gene.col] %in% genes.classed[[base.class]]), vec.space.axes])
            res <- cbind(data.frame(Family = family.name, mean.base.tiss.vers = (1 - 
                cosDiag(base.class.cloud$stat.vec)/sqrt(2)), stringsAsFactors = FALSE), 
                Reduce(cbind, lapply(names(g.c.i[evolved.classes]), function(evol.class) {
                  # Check whether evol.class has expression values for its genes:
                  if (g.c.i[[evol.class]]) {
                    e.g <- genes.classed[[evol.class]]
                    e.g.cloud <- statVectorCloud(expr.vecs[which(expr.vecs[, 
                      gene.col] %in% e.g), vec.space.axes])
                    evol.df <- data.frame(V1 = (1 - cosDiag(e.g.cloud$stat.vec)/sqrt(2)), 
                      V2 = rad2deg(acos(cosAngleVec(base.class.cloud$orth.on.diag.2.stat.vec, 
                        e.g.cloud$orth.on.diag.2.stat.vec))), V3 = distVectorClouds(base.class.cloud, 
                        e.g.cloud), stringsAsFactors = FALSE)
                    colnames(evol.df) <- paste(evol.class, c("tiss.vers", 
                      "tiss.change", "dist.vec.clouds"), sep = ".")
                    evol.df
                  } else {
                    evol.df <- data.frame(V1 = NA, V2 = NA, V3 = NA, stringsAsFactors = FALSE)
                    colnames(evol.df) <- paste(evol.class, c("tiss.vers", 
                      "tiss.change", "dist.vec.clouds"), sep = ".")
                    evol.df
                  }
                })))
        }
    }
    res
}


tands.w.orths


removeExpressionVariant <- function(gene.ids, reg.ex = getOption("GeneFamilies.expression.variant.regex", "\\.\\d+$")) {unique(sub(reg.ex, "", gene.ids))}


t.i <- which(as.logical(lapply(paralogs.lst, function(t.c) any(removeExpressionVariant(t.c) %in% orthologs.genes))))
para.w.orths <- lapply(paralogs.lst[t.i], function(t.c) removeExpressionVariant(t.c))


#' Investigate what happened in each gene group after duplication:
#' - Comparing tandems and orthologs
tands.w.orths.angles.df <- Reduce(rbind, mclapply(names(tands.w.orths), 
    function(fam.name) {
        genes <- intersect(tands.w.orths[[fam.name]], rpkm.expr.profiles.df$Gene)
        exprVecSpaceEvolutionAfterDupl(genes, fam.name, tand.classifier)
    }))


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
pdf(file.path(input.args[[1]], "inst", "meanTissueVersatilityDiffsAfterDuplicationBoxplot.pdf"))

par(old.par)

#' Plot both changes in tissue specificity:
p.df.1 <- dupl.w.orths.angles.df[with(dupl.w.orths.angles.df, which(!is.na(duplicated.tiss.change) & 
    !is.nan(duplicated.tiss.change))), ]
p.df.1$class <- "Duplicated"
p.df.2 <- tands.w.orths.angles.df[with(tands.w.orths.angles.df, which(!is.na(tandem.tiss.change) & 
    !is.nan(tandem.tiss.change))), ]
colnames(p.df.2) <- sub("tandem", "duplicated", colnames(p.df.2))
p.df.2$class <- "Tandem"
cols.i <- c("class")
p.df <- rbind(p.df.1, p.df.2)

#' Plot all
pdf(file.path(input.args[[1]], "inst", "meanChangeInAngleBetweenOrthOnDiagsBoxplot.pdf"))
pushed.colors <- brewer.pal(3, "Dark2")
boxplot(duplicated.tiss.change ~ class, data = p.df, xlab = "Type of Gene", 
    ylab = "angle between diagonal orthologs", border = pushed.colors, 
    col = addAlpha(pushed.colors))
dev.off()













#' Boxplot log2 fold change of ortholog, tandem, and duplicated expressed
#' genes. Gene Expression is inferred on the intra-species level during fruit
#' development and shade reaction.
