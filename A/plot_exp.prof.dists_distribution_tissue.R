
library(dotenv)
# Define directories for output data and results using environment variables
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Create results directory if it doesn't exist
if(!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}

message("USAGE: Rscript exec/plot_exp.prof.dists_distribution_tissue.R")

# Librarys for handling Dataframes, Lists etc. more efficiently
library(dplyr)
library(tidyr)
library(purrr)
library(parallel)

# Librarys necessary for plotting
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(ggpubr)

# ------------------------------------------------------------------------

# Automatically select data objects that contain tissue-specific statistics
load(file.path(output_data_dir, "exp.prof.dists_statistics.RData"))
loaded_objects <- ls()
data_names_tissue <- loaded_objects[grepl("_tissue_stats$", loaded_objects)]

# Process each tissue dataset to extract and organize mean and median distances by tissue and type
# Initialize empty data frames to hold tissue-specific mean and median distance data
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

# Generate boxplots for mean expression distances by tissue type
output_pdf <- file.path(results_dir, "tissues_mean_boxplot_combined_(t.test).pdf")
tissue_types <- unique(df_mean.dists_tissue$Tissue)

# Create empty list to store plots
plot_list <- list()

# Generate plot for each tissue
for (tissue in tissue_types) {
  df_tissue <- subset(df_mean.dists_tissue, Tissue == tissue)
 
  boxplot_tissue <- ggplot(df_tissue, aes(x = Type, y = Distance, fill = Type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    labs(title = paste("Mean Expression Distances (t.test) - ", tissue), y = "Distance", x = "Gene Type") +
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
      geom_signif(
        comparisons = type_combinations,
        test = "t.test",
        test.args = list(alternative = "two.sided"),
        map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
        step_increase = 0.05,
        tip_length = 0.005,
        color = "black",
        size = 0.3,
        textsize = 2.5
      )
  }
  
  # Add plot to list
  plot_list[[tissue]] <- boxplot_tissue
}

# Save all plots with minimal page numbers
ggsave(output_pdf, marrangeGrob(plot_list, nrow=1, ncol=1, top=""), width = 12, height = 8)

# ---------------------------------------------------------------------------

# Generate boxplots for median expression distances by tissue type
output_pdf <- file.path(results_dir, "tissues_median_boxplot_combined_(t.test).pdf")
tissue_types <- unique(df_median.dists_tissue$Tissue)

# Create empty list to store plots
plot_list <- list()

# Generate plot for each tissue
for (tissue in tissue_types) {
  df_tissue <- subset(df_median.dists_tissue, Tissue == tissue)
 
  boxplot_tissue <- ggplot(df_tissue, aes(x = Type, y = Distance, fill = Type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    labs(title = paste("Median Expression Distances (t.test) - ", tissue), y = "Distance", x = "Gene Type") +
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
      geom_signif(
        comparisons = type_combinations,
        test = "t.test",
        test.args = list(alternative = "two.sided"),
        map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
        step_increase = 0.05,
        tip_length = 0.005,
        color = "black",
        size = 0.3,
        textsize = 2.5
      )
  }
  
  # Add plot to list
  plot_list[[tissue]] <- boxplot_tissue
}

# Save all plots with minimal page numbers
ggsave(output_pdf, marrangeGrob(plot_list, nrow=1, ncol=1, top=""), width = 12, height = 8)


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------


# Generate boxplots for mean expression distances by tissue type
output_pdf <- file.path(results_dir, "tissues_mean_boxplot_combined_(wilcox.test).pdf")
tissue_types <- unique(df_mean.dists_tissue$Tissue)

# Create empty list to store plots
plot_list <- list()

# Generate plot for each tissue
for (tissue in tissue_types) {
  df_tissue <- subset(df_mean.dists_tissue, Tissue == tissue)
 
  boxplot_tissue <- ggplot(df_tissue, aes(x = Type, y = Distance, fill = Type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    labs(title = paste("Mean Expression Distances (wilcox.test) - ", tissue), y = "Distance", x = "Gene Type") +
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
      geom_signif(
        comparisons = type_combinations,
        test = "wilcox.test",
        test.args = list(alternative = "two.sided"),
        map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
        step_increase = 0.05,
        tip_length = 0.005,
        color = "black",
        size = 0.3,
        textsize = 2.5
      )
  }
  
  # Add plot to list
  plot_list[[tissue]] <- boxplot_tissue
}

# Save all plots with minimal page numbers
ggsave(output_pdf, marrangeGrob(plot_list, nrow=1, ncol=1, top=""), width = 12, height = 8)

# ---------------------------------------------------------------------------

# Generate boxplots for median expression distances by tissue type
output_pdf <- file.path(results_dir, "tissues_median_boxplot_combined_(wilcox.test).pdf")
tissue_types <- unique(df_median.dists_tissue$Tissue)

# Create empty list to store plots
plot_list <- list()

# Generate plot for each tissue
for (tissue in tissue_types) {
  df_tissue <- subset(df_median.dists_tissue, Tissue == tissue)
 
  boxplot_tissue <- ggplot(df_tissue, aes(x = Type, y = Distance, fill = Type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    labs(title = paste("Median Expression Distances (wilcox.test) - ", tissue), y = "Distance", x = "Gene Type") +
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
      geom_signif(
        comparisons = type_combinations,
        test = "wilcox.test",
        test.args = list(alternative = "two.sided"),
        map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
        step_increase = 0.05,
        tip_length = 0.005,
        color = "black",
        size = 0.3,
        textsize = 2.5
      )
  }
  
  # Add plot to list
  plot_list[[tissue]] <- boxplot_tissue
}

# Save all plots with minimal page numbers
ggsave(output_pdf, marrangeGrob(plot_list, nrow=1, ncol=1, top=""), width = 12, height = 8)

