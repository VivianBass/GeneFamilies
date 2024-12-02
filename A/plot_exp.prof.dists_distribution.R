
# Load environment variables to define directories for output data and results
library(dotenv)
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Create results directory if it doesn't exist
if(!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}

message("USAGE: Rscript exec/plot_exp.prof.dists_distribution.R")

# Libraries for efficient data handling
library(dplyr)
library(tidyr)
library(purrr)
library(parallel)

# Libraries for data visualization
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(ggpubr)

# ------------------------------------------------------------------------

# Load expression profile distance statistics, for each gene group (5 in total)
load(file.path(output_data_dir, "exp.prof.dists_statistics.RData"))
loaded_objects <- ls()
data_names <- loaded_objects[grepl(".lst_dists_stats$", loaded_objects)]

# Initialize empty data frames to store the mean and median distances for each dataset
# Header Type and Distance, in long format, Type containing the gene-groups
df_mean.dists <- data.frame()
df_median.dists <- data.frame()

# Process each dataset to extract and compile mean and median distance information
for (data_name in data_names) {
    current_data <- get(data_name)
    type_name <- sub(".lst_dists_stats", "", data_name)
    
    # Compile mean distance data into a tibble, adding columns for type, cluster, and distance
    temp_mean_df <- tibble(
        Type = type_name,
        Cluster = names(current_data$Mean),
        Distance = unlist(current_data$Mean)
    )
    df_mean.dists <- bind_rows(df_mean.dists, temp_mean_df)
    
    # Compile median distance data into a tibble, adding columns for type, cluster, and distance
    temp_median_df <- tibble(
        Type = type_name,
        Cluster = names(current_data$Median),
        Distance = unlist(current_data$Median)
    )
    df_median.dists <- bind_rows(df_median.dists, temp_median_df)
}

# Save both dataframes to the output directory
save(df_mean.dists, df_median.dists, 
     file = file.path(output_data_dir, "exp.prof.dists_mean_median.RData"))

# ---------------------------------------------------------------------------

# Define combinations of types for pairwise significance annotations in boxplots
types <- unique(df_mean.dists$Type)
type_combinations <- combn(types, 2, simplify = FALSE)


# Plot the boxplot with significance annotations for mean distances
boxplot_mean <- ggplot(df_mean.dists, aes(x = Type, y = Distance, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Mean Expression Distances (t.test)", y = "Distance") +
  theme_pubr(border = TRUE) +
  scale_y_continuous(breaks = seq(0, max(df_mean.dists$Distance), by = 0.2)) +
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
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
    axis.text.x = element_text(size = 10),
    plot.margin = margin(r = 30)
  )

boxplot_median <- ggplot(df_median.dists, aes(x = Type, y = Distance, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Median Expression Distances (t.test)", y = "Distance") +
  theme_pubr(border = TRUE) +
  scale_y_continuous(breaks = seq(0, max(df_median.dists$Distance), by = 0.2)) +
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
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
    axis.text.x = element_text(size = 10),   
    plot.margin = margin(r = 30)
  )
# Save both plots to a multi-page PDF
output_pdf <- file.path(results_dir, "boxplots_expression_distances_(t.test).pdf")

# Combine plots into a list
plot_list_tea <- list(boxplot_median, boxplot_mean)

# Save all plots with minimal page numbers
ggsave(output_pdf, marrangeGrob(plot_list_tea, nrow=1, ncol=1, top=""), 
       width = 12, height = 8, device = "pdf")


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------


# Define combinations of types for pairwise significance annotations in boxplots
types <- unique(df_mean.dists$Type)
type_combinations <- combn(types, 2, simplify = FALSE)


# Plot the boxplot with significance annotations for mean distances
boxplot_mean <- ggplot(df_mean.dists, aes(x = Type, y = Distance, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Mean Expression Distances (wilcox.test)", y = "Distance") +
  theme_pubr(border = TRUE) +
  scale_y_continuous(breaks = seq(0, max(df_mean.dists$Distance), by = 0.2)) +
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
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
    axis.text.x = element_text(size = 10),
    plot.margin = margin(r = 30)
  )

boxplot_median <- ggplot(df_median.dists, aes(x = Type, y = Distance, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Median Expression Distances (wilcox.test)", y = "Distance") +
  theme_pubr(border = TRUE) +
  scale_y_continuous(breaks = seq(0, max(df_median.dists$Distance), by = 0.2)) +
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
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
    axis.text.x = element_text(size = 10),   
    plot.margin = margin(r = 30)
  )
# Save both plots to a multi-page PDF
output_pdf <- file.path(results_dir, "boxplots_expression_distances_(wilcox.test).pdf")

# Combine plots into a list
plot_list_wilcox <- list(boxplot_median, boxplot_mean)

# Save all plots with minimal page numbers
ggsave(output_pdf, marrangeGrob(plot_list_wilcox, nrow=1, ncol=1, top=""), 
       width = 12, height = 8, device = "pdf")