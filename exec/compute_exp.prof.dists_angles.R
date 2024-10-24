
require(GeneFamilies)
options(mc.cores = getMcCores())

library(parallel)
library(dotenv)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

source("R/expression_funks.R")
load(file.path(output_data_dir, "RNA_Seq_RPKM_and_profiles.RData"))

# Define expression vector space (you can modify this as needed)
# The number of tissues corresponds to the count of axes, excluding the gene column in the names.
expr.cols <- names(rpkm.expr.profiles.df)[names(rpkm.expr.profiles.df) != "Gene"]
n.dims <- length(expr.cols)

message("USAGE: Rscript exec/compute_exp.prof.dists_angles.R")
message("Defining expression vector space with the following axes: ", paste(expr.cols, collapse = ", "))

source("R/expression_funks.R")
load(file.path(output_data_dir, "RNA_Seq_RPKM_and_profiles.RData"))


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


# merge paralog.expr.angle.diag.df and orths.expr.angle.diag.df and plot results:

p.lst <- list(paralog = paralog.expr.angle.diag.df, ortholog = orths.expr.angle.diag.df)
expr.angle.diag.df <- Reduce(rbind, mclapply(names(p.lst), function(gene.type) {
    data.frame(gene.type = gene.type, angle.diag = p.lst[[gene.type]]$angle.diag, 
        stringsAsFactors = FALSE)
}))

save(p.lst, para.expr.angle.diag.df, orths.expr.angle.diag.df, expr.angle.diag.df, 
    file = file.path(output_data_dir, "expr.angle.diag.RData"))

message("DONE")











