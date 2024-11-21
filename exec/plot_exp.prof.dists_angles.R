
require(GeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript exec/plot_exp.prof.dists_angles.R")

library(parallel)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(dotenv)

# Set-up output directory, defined in the .env file 
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# load gene-groups angles datasets
load(file.path(output_data_dir, "exp.prof.dists_angles.RData"))

# functions sourced from:
source("R/angles_funks.R")


# Create list of dataframes to validate
df_list <- list(
    con_orthologs = con_orthologs.expr.angle.diag.df,
    in_paralogs = in_paralogs.expr.angle.diag.df,
    out_paralogs = out_paralogs.expr.angle.diag.df,
    special_in_paralogs = special_in_paralogs.expr.angle.diag.df,
    special_out_paralogs = special_out_paralogs.expr.angle.diag.df
)

# Create p.lst with only valid dataframes while preserving names
p.lst <- validate_angle_dataframes(df_list)

p.df <- Reduce(rbind, mclapply(names(p.lst), function(gene.type) {
    data.frame(
        gene.type = gene.type,
        angle.diag = p.lst[[gene.type]]$angle.diag,
        stringsAsFactors = FALSE
    )
}))

# Get actual unique values from data
actual_levels <- unique(p.df$gene.type)

plot.df <- p.df[!is.nan(p.df$angle.diag), ]
plot.df$gene.type <- factor(plot.df$gene.type, levels = actual_levels)

pdf(file.path(results_dir, "expressionAngleToDiagonalBoxplot.pdf"))
colors <- brewer.pal(length(p.lst), "Dark2")
pushed.colors <- append(colors[3], colors[1:2])
boxplot(angle.diag ~ gene.type, data = plot.df, 
        xlab = "Type of Gene", 
        ylab = "relative tissue specificity", 
        border = pushed.colors, 
        col = addAlpha(pushed.colors))
dev.off()


# Calculate relative versatility
plot.df$rel.vers <- 1 - plot.df$angle.diag

# Use actual gene types from data
actual_levels <- unique(plot.df$gene.type)
plot.df$gene.type <- factor(plot.df$gene.type, levels = actual_levels)

pdf(file.path(results_dir, "relativeExpressionVersatilityBoxplot.pdf"))
colors <- brewer.pal(length(p.lst), "Dark2")
pushed.colors <- append(colors[3], colors[1:2])
boxplot(rel.vers ~ gene.type, data = plot.df, 
        xlab = "Type of Gene", 
        ylab = "relative tissue versatility",
        border = pushed.colors, 
        col = addAlpha(pushed.colors), 
        outline = FALSE)
dev.off()



