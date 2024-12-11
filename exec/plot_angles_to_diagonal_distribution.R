
# Import necessary libraries
library(dotenv) # Load environment variables for directory paths
library(dplyr) # Data manipulation
library(tidyr) # Data tidying
library(purrr) # Functional programming for data
library(tibble) # Handling tidy data frames
library(parallel) # Parallel processing for performance optimization
library(ggplot2) # Data visualization
library(ggsignif) # Add significance annotations to ggplots
library(gridExtra) # Arrange multiple plots
library(ggpubr) # Publication-ready ggplots
library(rstatix) # Perform statistical tests (t-tests, Wilcoxon tests)

# Load user-defined functions
source("R/compute_funks.R") # Functions for computing distance statistics
source("R/angles_to_diagonal_funks.R")

# Set output directories from environment variables
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Ensure the results directory exists, create it if not
if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}

message("USAGE: Rscript exec/plot_angles_to_diagonal_distribution.R")

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

# Call the function with the list of test types
create_angle_versatility_plots(plot.df, c("t.test", "wilcox.test"), results_dir)


# ------------------------------------------------------------------------

# Load gene-groups angles datasets
load(file.path(output_data_dir, "exp.prof.dists_angles_log2.RData"))

# Create a list of dataframes to validate
df_list_log2 <- list(
    con_orthologs = con_orthologs.expr.angle.diag_log2.df,
    in_paralogs = in_paralogs.expr.angle.diag_log2.df,
    out_paralogs = out_paralogs.expr.angle.diag_log2.df,
    special_in_paralogs = special_in_paralogs.expr.angle.diag_log2.df,
    special_out_paralogs = special_out_paralogs.expr.angle.diag_log2.df
)

# Validate and filter dataframes
p.lst_log2 <- validate_angle_dataframes(df_list_log2)

# Combine validated data into a single dataframe
p.df_log2 <- Reduce(rbind, mclapply(names(p.lst_log2), function(gene.type) {
    data.frame(
        gene.type = gene.type,
        angle.diag = p.lst_log2[[gene.type]]$angle.diag,
        stringsAsFactors = FALSE
    )
}))

# Filter and prepare data
plot.df_log2 <- p.df_log2[!is.nan(p.df_log2$angle.diag), ]
plot.df_log2$gene.type <- factor(plot.df_log2$gene.type, levels = unique(plot.df_log2$gene.type))
plot.df_log2$rel.vers <- 1 - plot.df_log2$angle.diag

# Setup
gene_types <- levels(plot.df_log2$gene.type)
type_combinations <- combn(gene_types, 2, simplify = FALSE)

# Call the function with the list of test types and specify log2 in the output
create_angle_versatility_plots(plot.df_log2, c("t.test", "wilcox.test"), results_dir, suffix = "_log2")

# --------------------------------------------------------------------------------

# Perform statistical tests for columns "angle.diag" "rel.vers" in both dataframes
test_results <- perform_statistical_tests_for_columns(
    plot.df,
    plot.df_log2,
    "statistical_tests_angle_to_diagonal_and_rel.vers.csv"
)




