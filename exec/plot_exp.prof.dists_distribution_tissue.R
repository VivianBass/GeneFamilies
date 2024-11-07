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

significance_level <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")  
}

# Automatically sort the loadedstatistics data into regular and tissue datasets
load(file.path(output_data_dir, "exp.prof.dists_statistics.RData"))
loaded_objects <- ls()
data_names_tissue <- loaded_objects[grepl("\\.dists\\.tissue_stats$", loaded_objects)]

# Initialize empty dataframes for tissue data
df_mean.dists_tissue <- data.frame()
df_median.dists_tissue <- data.frame()

# Process each tissue dataset
for (data_name in data_names_tissue) {
    current_data <- get(data_name)
    
    type_name <- sub("_filtered.filtered.dists.tissue_stats", "", data_name)
    
    # Create and combine mean data
    temp_mean_df <- tibble(
        Type = type_name,
        Tissue = names(current_data$Mean),
        Distance = unlist(current_data$Mean)
    )
    df_mean.dists_tissue <- bind_rows(df_mean.dists_tissue, temp_mean_df)
    
    # Create and combine median data
    temp_median_df <- tibble(
        Type = type_name,
        Tissue = names(current_data$Median),
        Distance = unlist(current_data$Median)
    )
    df_median.dists_tissue <- bind_rows(df_median.dists_tissue, temp_median_df)
}

# Perform t-tests for median and mean per tissue between orthologs and paralogs group 
# and adjust p-values, then apply the significance_level function
# group_by() -> defines the tissue column

t_test_median_tissue <- df_median.dists_tissue %>%
  group_by(variables) %>%                             
  t_test(value ~ Type, alternative = "greater") %>%  
  adjust_pvalue(method = "BH") %>%                   
  mutate(significance = sapply(p, significance_level),  
         analysis = "Median")  

t_test_mean_tissue <- df_mean.dists_tissue %>%
  group_by(variables) %>%                             
  t_test(value ~ Type, alternative = "greater") %>%  
  adjust_pvalue(method = "BH") %>%                   
  mutate(significance = sapply(p, significance_level),  
         analysis = "Mean")  

t_test_summary_tissue <- bind_rows(t_test_median_tissue, t_test_mean_tissue)

write.csv(t_test_summary_tissue, file.path(results_dir, "t_test_summary_tissue.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------

# geom_signif(): This function is used to add significance annotations directly to the plot 
# so plots already have the significance level !!

output_pdf <- file.path(results_dir, "tissues_mean_boxplot_combined.pdf")
# Get unique tissue types from the dataset
tissue_types <- unique(df_mean.dists_tissue$Tissue)

# Open the PDF device
pdf(output_pdf, width = 12, height = 8)

# Loop through each tissue type and create a plot
for (tissue in tissue_types) {
  # Filter data for the current tissue type
  df_tissue <- subset(df_mean.dists_tissue, Tissue == tissue)
  
  # Create the boxplot for the current tissue type
  boxplot_tissue <- ggplot(df_tissue, aes(x = Type, y = Distance, fill = Type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    labs(title = paste("Mean Expression Distances -", tissue), y = "Distance", x = "Gene Type") +
    theme_pubr(border = TRUE) +
    theme(
      plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
      axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
      axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20))
    )
  
  # Only add statistical significance if both comparison groups are present
  if (all(c("ortholog", "paralog") %in% df_tissue$Type)) {
    boxplot_tissue <- boxplot_tissue +
      geom_signif(comparisons = list(c("ortholog", "paralog")), map_signif_level = TRUE)
  }
  # Print the plot to the PDF file
  print(boxplot_tissue)
} 
dev.off()

df_median.dists_tissue
# ---------------------------------------------------------------------------

# Creating a Boxplot for Median Distances
# Calculate the appropriate height for the PDF based on the number of facets
output_pdf <- file.path(results_dir, "tissues_median_boxplot_combined.pdf")
# Get unique tissue types from the dataset
tissue_types <- unique(df_median.dists_tissue$Tissue)

# Open the PDF device
pdf(output_pdf, width = 12, height = 8)

# Loop through each tissue type and create a plot
for (tissue in tissue_types) {
  # Filter data for the current tissue type
  df_tissue <- subset(df_median.dists_tissue, Tissue == tissue)
  
  # Create the boxplot for the current tissue type
  boxplot_tissue <- ggplot(df_tissue, aes(x = Type, y = Distance, fill = Type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    labs(title = paste("Median Expression Distances -", tissue), y = "Distance", x = "Gene Type") +
    theme_pubr(border = TRUE) +
    theme(
      plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
      axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
      axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20))
    )
  
  # Only add statistical significance if both comparison groups are present
  if (all(c("Ortholog", "Paralog") %in% df_tissue$Type)) {
    boxplot_tissue <- boxplot_tissue +
      geom_signif(comparisons = list(c("Ortholog", "Paralog")), map_signif_level = TRUE)
  }
  
  # Print the plot to the PDF file
  print(boxplot_tissue)
}
dev.off()


