
# Load environment variables to define directories for output data and results
library(dotenv)
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Create results directory if it doesn't exist
if(!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}

message("USAGE: Rscript exec/plot_exp.prof.dists_angles.R")

# Libraries for efficient data handling
library(dplyr)
library(tidyr)
library(purrr)
library(parallel)

# Libraries for data visualization
library(RColorBrewer)
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(ggpubr)

# Source custom functions
source("R/angles_funks.R")

# ------------------------------------------------------------------------

# Load gene-groups angles datasets
load(file.path(output_data_dir, "exp.prof.dists_angles.RData"))

# Create a list of dataframes to validate
df_list <- list(
    con_orthologs = con_orthologs.expr.angle.diag.df,
    in_paralogs = in_paralogs.expr.angle.diag.df,
    out_paralogs = out_paralogs.expr.angle.diag.df,
    special_in_paralogs = special_in_paralogs.expr.angle.diag.df,
    special_out_paralogs = special_out_paralogs.expr.angle.diag.df
)

# Validate and filter dataframes
p.lst <- validate_angle_dataframes(df_list)

# Combine validated data into a single dataframe
p.df <- Reduce(rbind, mclapply(names(p.lst), function(gene.type) {
    data.frame(
        gene.type = gene.type,
        angle.diag = p.lst[[gene.type]]$angle.diag,
        stringsAsFactors = FALSE
    )
}))

# Filter and prepare data
plot.df <- p.df[!is.nan(p.df$angle.diag), ]
plot.df$gene.type <- factor(plot.df$gene.type, levels = unique(plot.df$gene.type))
plot.df$rel.vers <- 1 - plot.df$angle.diag

# Setup
gene_types <- levels(plot.df$gene.type)
type_combinations <- combn(gene_types, 2, simplify = FALSE)

# Generate all plots
for (test in c("t.test", "wilcox.test")) {
    create_angle_versatility_plot(plot.df, "angle", test)
    create_angle_versatility_plot(plot.df, "versatility", test)
}