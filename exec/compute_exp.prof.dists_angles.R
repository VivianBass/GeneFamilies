
require(GeneFamilies)
options(mc.cores = getMcCores())

library(parallel)
library(dotenv)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")

message("USAGE: Rscript path/2/GeneFamilies/exec/expressionAngles.R path/2/GeneFamilies")
message("Defining expression vector space with the following axes: ", paste(expr.cols, collapse = ", "))

source("R/expression_funks.R")

load(file.path(output_data_dir, "RNA_Seq_RPKM_and_profiles.RData"))



# need the expression lists (median) for paralogs and orthologs



load("data/data.RData")

# Inputs used in this script:
rna.seq.exp.profils <- df14_P_dmel
para.expr.lst <- vv_df3_para.lst
orths.expr.lst <- vv_df3_ortho.lst


input.args <- list()
input.args[[1]] <- orths.expr.lst
input.args[[2]] <- para.expr.lst
input.args[[3]] <- rna.seq.exp.profils

# Define expression vector space (you can modify this as needed)
expr.cols <- c("Whole_Body", "Muscle", "Gut", "Fat_Body")
n.dims <- length(expr.cols)

# - need a function for this


# calculate angles between expression profiles and diagonal , as measure of tissue specificity
# Function to calculate angles to diagonal for a list of genes
calculate_angle_to_diagonal <- function(gene_list, expression_profiles, expression_columns) {
  # Find intersection of input gene list with available genes in the expression profile data
  genes_with_expr <- intersect(unlist(gene_list), expression_profiles$Gene)
  
  # Calculate the angle to the diagonal for each gene
  angle_diag_df <- data.frame(
    Gene = genes_with_expr,
    angle.diag = as.numeric(mclapply(genes_with_expr, function(gene) {
      # Compute cosine diagonal distance for the expression profile of the gene
      cosDiag(expression_profiles[which(expression_profiles$Gene == gene), expression_columns]) / sqrt(2)
    })),
    stringsAsFactors = FALSE
  )
  
  return(angle_diag_df)
}

para.expr.angle.diag.df <- calculate_angle_to_diagonal(input.args[[1]], input.args[[3]], expr.cols)
orths.expr.angle.diag.df <- calculate_angle_to_diagonal(input.args[[2]], input.args[[3]], expr.cols)

# merge paralog.expr.angle.diag.df and orths.expr.angle.diag.df as one df:
p.lst <- list(paralog = paralog.expr.angle.diag.df, ortholog = orths.expr.angle.diag.df)
expr.angle.diag.df <- Reduce(rbind, mclapply(names(p.lst), function(gene.type) {
    data.frame(gene.type = gene.type, angle.diag = p.lst[[gene.type]]$angle.diag, 
        stringsAsFactors = FALSE)
}))

save(para.expr.angle.diag.df, orths.expr.angle.diag.df, expr.angle.diag.df,
    file = file.path(output_data_dir, "expr.angle.diag.RData"))

message("DONE")

# - would need to clear for invalid numbers









