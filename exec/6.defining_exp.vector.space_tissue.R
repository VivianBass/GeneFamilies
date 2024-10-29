
require(GeneFamilies)
options(mc.cores = getMcCores())
library(dotenv)

library(ggplot2)
library(ggsignif)
library(gridExtra)

library(dplyr)
library(purrr)

library(tidyverse)
library(rstatix)
library(ggpubr)
library(tibble)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

load(file.path(output_data_dir, "exp.prof.dists.RData"))


# wich dataframes do we start of here ???

# decide naming convention of the dataframes 

# Define expression vector space:
# get the tissues directly from the expression profiles 

# Filter the data to include only the selected tissues, should this step be done before computation or after computation ???
# Genegroups used for the analysis are: Orthologs and Paralogs
# Filtering and preprocessing the data for the required statistics: mean and median


selected_tissues <- c("mE_mRNA_1182-4H_cells", "mE_mRNA_A_1d_carcass", "mE_mRNA_A_1d_dig_sys","mE_mRNA_A_20d_carcass","mE_mRNA_A_20d_dig_sys") 

filtered_data <- boxplot_data %>% filter(tissue %in% selected_tissues)
save(filtered_data, file = file.path(output_data_dir, "/filtered_data_expression_distances.RData"))
filtered_data <- na.omit(filtered_data)

# --------------------------------------------------------------------------

# assuming already have the selected tissues and the mean / median expre values
# get the colnames 

expr.cols.para <- names(df_median_mean_paralogs.exp.prof.dists.tissue)
expr.cols.ortho <- names(df_median_mean_orthologs.exp.prof.dists.tissue)

expr.cols <- c("Whole_Body", "Muscle", "Gut", "Fat_Body")
n.dims <- length(expr.cols)

# --------------------------------------------------------------------------

# Create a Combined Mean and Median DataFrame for the Orthologs and Paralogs Gene Groups


df_ortho_mean_tissue <- df_median_mean_paralogs.exp.prof.dists.tissue[ ,  c("Cluster", "Mean_Fat_Body", "Mean_Gut", "Mean_Muscle", "Mean_Whole_Body")]

df_para_mean_tissue <- df_median_mean_paralogs.exp.prof.dists.tissue[ ,  c("Cluster", "Mean_Fat_Body", "Mean_Gut", "Mean_Muscle", "Mean_Whole_Body")]

df_ortho_median_tissue <- df_median_mean_orthologs.exp.prof.dists.tissue[ ,  c("Cluster", "Median_Fat_Body", "Median_Gut", "Median_Muscle", "Median_Whole_Body")]

df_para_median_tissue <- df_median_mean_orthologs.exp.prof.dists.tissue[ ,  c("Cluster", "Median_Fat_Body", "Median_Gut", "Median_Muscle", "Median_Whole_Body")]


mean_p.lst <- list(Paralog = df_ortho_mean_tissue, Ortholog = df_para_mean_tissue)
mean_col_names <- names(mean_p.lst[[1]])[names(mean_p.lst[[1]]) != "Cluster"]
df_ortho_para_mean_tissue <- map_dfr(mean_p.lst, ~as.data.frame(.x[mean_col_names]), .id = "Type")

median_p.lst <- list(Paralog = df_para_median_tissue, Ortholog = df_ortho_median_tissue)
median_col_names <- names(median_p.lst[[1]])[names(median_p.lst[[1]]) != "Cluster"]
df_ortho_para_median_tissue <- map_dfr(median_p.lst, ~as.data.frame(.x[median_col_names]), .id = "Type")


# Converting DataFrames to long format is ideal for plotting, 
# as it helps avoid issues with setting the size of DataFrames if datasets differ in length.

mean_p.df_long <- mean_p.df %>% pivot_longer(cols = -Type, names_to = "variables", values_to = "value") %>% 
                  filter(!is.na(value) & !is.infinite(value))

median_p.df_long <- median_p.df %>% pivot_longer(cols = -Type, names_to = "variables", values_to = "value") %>% 
                  filter(!is.na(value) & !is.infinite(value))

save(df_ortho_para_mean_tissue, df_ortho_para_median_tissue, mean_p.df_long, median_p.df_long, 
     file = file.path(output_data_dir, "exp.prof.dists_tissue.RData"))