require(GeneFamilies)
options(mc.cores = getMcCores())
library(parallel)

message("USAGE: Rscript exec/plot_exp.prof.dists_distribution_tissue.R")

library(dotenv)
# Define directories for output data and results using environment variables
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Librarys for handling Dataframes, Lists etc. more efficiently
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

# Librarys necessary for plotting
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(rstatix)
library(ggpubr)

# ------------------------------------------------------------------------

# Automatically select data objects that contain tissue-specific statistics
load(file.path(output_data_dir, "exp.prof.dists_statistics.RData"))
loaded_objects <- ls()
data_names_tissue <- loaded_objects[grepl("_tissue_stats$", loaded_objects)]

# Initialize empty data frames to hold tissue-specific mean and median distance data
df_mean.dists_tissue <- data.frame()
df_median.dists_tissue <- data.frame()

# Process each tissue dataset to extract and organize mean and median distances by tissue and type
# Initialize empty data frames
df_mean.dists_tissue <- data.frame()
df_median.dists_tissue <- data.frame()

# Process each tissue dataset
for (data_name in data_names_tissue) {
    current_data <- get(data_name)  
    type_name <- sub(".lst_dists_tissue_stats", "", data_name)  
    
    # Extract data and arrange columns with Type in second position
    temp_mean_df <- current_data %>%
        select(Family, Tissue, Mean) %>%
        rename(Distance = Mean) %>%
        mutate(Type = type_name) %>%
        select(Family, Type, Tissue, Distance)
    
    temp_median_df <- current_data %>%
        select(Family, Tissue, Median) %>%
        rename(Distance = Median) %>%
        mutate(Type = type_name) %>%
        select(Family, Type, Tissue, Distance)
    
    df_mean.dists_tissue <- bind_rows(df_mean.dists_tissue, temp_mean_df)
    df_median.dists_tissue <- bind_rows(df_median.dists_tissue, temp_median_df)
}

# Save both dataframes to the output directory
save(df_mean.dists_tissue, df_median.dists_tissue, 
     file = file.path(output_data_dir, "exp.prof.dists_mean_median_tissue.RData"))

# ---------------------------------------------------------------------------

# Function to assign significance level based on p-value
significance_level <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")  
}

# ---------------------------------------------------------------------------
df_mean.dists_tissue
# Generate boxplots for mean expression distances by tissue type
output_pdf <- file.path(results_dir, "tissues_mean_boxplot_combined.pdf")
tissue_types <- unique(df_mean.dists_tissue$Tissue)

pdf(output_pdf, width = 12, height = 8)

for (tissue in tissue_types) {
  df_tissue <- subset(df_mean.dists_tissue, Tissue == tissue)
  
  boxplot_tissue <- ggplot(df_tissue, aes(x = Type, y = Distance, fill = Type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    labs(title = paste("Mean Expression Distances -", tissue), y = "Distance", x = "Gene Type") +
    theme_pubr(border = TRUE) +
    scale_y_continuous(breaks = seq(0, max(df_tissue$Distance), by = 0.1)) +
    theme(
      plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
      axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
      axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
      axis.text.x = element_text(size = 10),  
      plot.margin = margin(r = 30)
    )
  
  # Add significance annotations only if there are multiple types
  types <- unique(df_tissue$Type)
  if (length(types) >= 2) {
    type_combinations <- combn(types, 2, simplify = FALSE)
    boxplot_tissue <- boxplot_tissue +
      geom_signif(comparisons = type_combinations, map_signif_level = TRUE)
  }
  
  print(boxplot_tissue)  # Output the plot to the PDF
}

dev.off()


# ---------------------------------------------------------------------------

# Generate boxplots for median expression distances by tissue type
output_pdf <- file.path(results_dir, "tissues_median_boxplot_combined.pdf")
tissue_types <- unique(df_median.dists_tissue$Tissue)
pdf(output_pdf, width = 12, height = 8)

for (tissue in tissue_types) {
  df_tissue <- subset(df_median.dists_tissue, Tissue == tissue)
  
  boxplot_tissue <- ggplot(df_tissue, aes(x = Type, y = Distance, fill = Type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    labs(title = paste("Median Expression Distances -", tissue), y = "Distance", x = "Gene Type") +
    theme_pubr(border = TRUE) +
    scale_y_continuous(breaks = seq(0, max(df_tissue$Distance), by = 0.1)) +
    theme(
      plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
      axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
      axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
      axis.text.x = element_text(size = 10),  
      plot.margin = margin(r = 30)
    )
  
  # Add significance annotations only if there are multiple types
  types <- unique(df_tissue$Type)
  if (length(types) >= 2) {
    type_combinations <- combn(types, 2, simplify = FALSE)
    boxplot_tissue <- boxplot_tissue +
      geom_signif(comparisons = type_combinations, map_signif_level = TRUE)
  }
  
  print(boxplot_tissue)  # Output the plot to the PDF
}
dev.off()
