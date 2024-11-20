


#' Calculate expression profiles distances

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


