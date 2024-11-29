
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
library(tibble)
library(parallel)

# Libraries for data visualization
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(ggpubr)

# ---------------------------------------------------------------------------

# Automatically select data objects
load(file.path(output_data_dir, "exp.prof.dists.RData"))
loaded_objects <- ls()
data_names <- loaded_objects[grepl(".lst_dists$", loaded_objects)]
data_names_tissue <- loaded_objects[grepl(".lst_dists_tissue$", loaded_objects)]

# Process regular distances
df_complete_dists <- data_names %>%
    map_df(function(data_name) {
        current_data <- get(data_name)
        type_name <- sub(".lst_dists", "", data_name)
        
        tibble(
            Type = type_name,
            Distance = unlist(current_data)
        )
    })

# Process tissue-specific distances
df_complete_dists_tissue <- data_names_tissue %>%
    map_df(function(data_name) {
        current_data <- get(data_name)
        type_name <- sub(".lst_dists_tissue", "", data_name)
        
        enframe(current_data, name = "Family") %>%
        unnest_longer(value) %>%
        unnest_longer(value) %>%
        rename(
            Tissue = value_id,
            Distance = value
        ) %>%
        mutate(Type = type_name) %>%
        select(Family, Type, Tissue, Distance)
    })

# ---------------------------------------------------------------------------

df_complete_dists_filtered <- df_complete_dists %>%
    filter(is.finite(Distance))

types <- unique(df_complete_dists_filtered$Type)
type_combinations <- combn(types, 2, simplify = FALSE)

# Calculate y-axis limits with filtered data
y_max <- max(df_complete_dists_filtered$Distance, na.rm = TRUE)
y_breaks <- seq(0, y_max, by = 0.2)

boxplot_complete_tea <- ggplot(df_complete_dists_filtered, aes(x = Type, y = Distance, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Complete Expression Distances (t.test)", y = "Distance") +
  theme_pubr(border = TRUE) +
  scale_y_continuous(breaks = y_breaks) +
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

boxplot_complete_wilcox <- ggplot(df_complete_dists_filtered, aes(x = Type, y = Distance, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Complete Expression Distances (wilcox.test)", y = "Distance") +
  theme_pubr(border = TRUE) +
  scale_y_continuous(breaks = y_breaks) +
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
output_pdf <- file.path(results_dir, "complete_boxplots_expression_distances_(t.test_wilcox).pdf")

# Combine plots into a list
plot_list <- list(boxplot_complete_tea, boxplot_complete_wilcox)

# Save all plots with minimal page numbers
ggsave(output_pdf, marrangeGrob(plot_list, nrow=1, ncol=1, top=""),
       width = 12, height = 8, device = "pdf")

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

# Filter non-finite values
df_complete_dists_tissue_filtered <- df_complete_dists_tissue %>%
    filter(is.finite(Distance))

# Get unique tissue types
tissue_types <- unique(df_complete_dists_tissue_filtered$Tissue)

# Create empty lists to store plots for both test types
plot_list_tea <- list()
plot_list_wilcox <- list()

# Generate plots for each tissue
for (tissue in tissue_types) {
    df_tissue <- subset(df_complete_dists_tissue_filtered, Tissue == tissue)
    
    # Calculate y-axis breaks for current tissue
    y_breaks <- seq(0, max(df_tissue$Distance, na.rm = TRUE), by = 0.2)
    
    # Base plot for t-test
    boxplot_tissue_tea <- ggplot(df_tissue, aes(x = Type, y = Distance, fill = Type)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
        labs(title = paste("Complete Expression Distances (t.test) -", tissue), 
             y = "Distance", x = "Gene Type") +
        theme_pubr(border = TRUE) +
        scale_y_continuous(breaks = y_breaks) +
        theme(
            plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
            axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
            axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
            axis.text.x = element_text(size = 10),
            plot.margin = margin(r = 30)
        )
    
    # Base plot for wilcoxon test
    boxplot_tissue_wilcox <- boxplot_tissue_tea +
        labs(title = paste("Complete Expression Distances (wilcox.test) -", tissue))
    
    # Add significance annotations if multiple types exist
    types <- unique(df_tissue$Type)
    if (length(types) >= 2) {
        type_combinations <- combn(types, 2, simplify = FALSE)
        
        # Add t-test annotations
        boxplot_tissue_tea <- boxplot_tissue_tea +
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
        
        # Add wilcoxon test annotations
        boxplot_tissue_wilcox <- boxplot_tissue_wilcox +
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
    
    # Add plots to respective lists
    plot_list_tea[[tissue]] <- boxplot_tissue_tea
    plot_list_wilcox[[tissue]] <- boxplot_tissue_wilcox
}

# Save t-test plots
output_pdf_tea <- file.path(results_dir, "complete_boxplots_tissues_expression_distances_(t.test).pdf")
ggsave(output_pdf_tea, marrangeGrob(plot_list_tea, nrow=1, ncol=1, top=""), 
       width = 12, height = 8)

# Save wilcoxon test plots
output_pdf_wilcox <- file.path(results_dir, "complete_boxplots_tissues_expression_distances_(wilcox.test).pdf")
ggsave(output_pdf_wilcox, marrangeGrob(plot_list_wilcox, nrow=1, ncol=1, top=""), 
       width = 12, height = 8)




