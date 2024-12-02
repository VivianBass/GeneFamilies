
# Load environment variables to define directories for output data and results
library(dotenv)
output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")

# Create results directory if it doesn't exist
if(!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}

message("USAGE: Rscript exec/plot_exp.prof.angles_distribution.R")

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

# ------------------------------------------------------------------------

# Automatically select data objects
load(file.path(output_data_dir, "exp.prof.angles.RData"))
loaded_objects <- ls()
data_names <- loaded_objects[grepl(".lst_cos_angles_dists$", loaded_objects)]

# Process cosine angles distances
df_complete_angles <- data_names %>%
    map_df(function(data_name) {
        current_data <- get(data_name)
        type_name <- sub(".lst_cos_angles_dists", "", data_name)
        
        tibble(
            Type = type_name,
            Angle = unlist(current_data)
        )
    })

# ---------------------------------------------------------------------------

# Filter non-finite values and calculate breaks
df_complete_angles_filtered <- df_complete_angles %>%
    filter(is.finite(Angle))

y_max <- max(df_complete_angles_filtered$Angle, na.rm = TRUE)
y_breaks <- seq(0, y_max, by = 0.2)

# Define combinations for significance testing
types <- unique(df_complete_angles_filtered$Type)
type_combinations <- combn(types, 2, simplify = FALSE)

# Plot with filtered data
boxplot_tea <- ggplot(df_complete_angles_filtered, aes(x = Type, y = Angle, fill = Type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    labs(title = "Complete Expression Angles (t.test)", y = "Angle") +
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

boxplot_wilcox <- ggplot(df_complete_angles_filtered, aes(x = Type, y = Angle, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  labs(title = "Complete Expression Angles (wilcox.test)", y = "Angle") +
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
output_pdf <- file.path(results_dir, "complete_boxplots_expression_angles_(t.test_wilcox).pdf")

# Combine plots into a list
plot_list <- list(boxplot_tea, boxplot_wilcox)

# Save all plots with minimal page numbers
ggsave(output_pdf, marrangeGrob(plot_list, nrow=1, ncol=1, top=""),
       width = 12, height = 8, device = "pdf")

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

# Automatically select data objects
load(file.path(output_data_dir, "exp.prof.angles.log2.RData"))
loaded_objects <- ls()
data_names <- loaded_objects[grepl(".lst_cos_angles_dists_log2$", loaded_objects)]

# Process cosine angles distances
df_complete_angles <- data_names %>%
    map_df(function(data_name) {
        current_data <- get(data_name)
        type_name <- sub(".lst_cos_angles_dists_log2", "", data_name)
        
        tibble(
            Type = type_name,
            Angle = unlist(current_data)
        )
    })

# ---------------------------------------------------------------------------

# Filter non-finite values and calculate breaks
df_complete_angles_filtered <- df_complete_angles %>%
    filter(is.finite(Angle))

y_max <- max(df_complete_angles_filtered$Angle, na.rm = TRUE)
y_breaks <- seq(0, y_max, by = 0.2)

# Define combinations for significance testing
types <- unique(df_complete_angles_filtered$Type)
type_combinations <- combn(types, 2, simplify = FALSE)

# Plot with filtered data - t-test
boxplot_tea <- ggplot(df_complete_angles_filtered, aes(x = Type, y = Angle, fill = Type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    labs(title = "Complete Expression Angles log2 (t.test)", y = "Angle") +
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

# Plot with filtered data - wilcoxon test
boxplot_wilcox <- ggplot(df_complete_angles_filtered, aes(x = Type, y = Angle, fill = Type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    labs(title = "Complete Expression Angles log2 (wilcox.test)", y = "Angle") +
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
output_pdf <- file.path(results_dir, "complete_boxplots_expression_angles_log2_(t.test_wilcox).pdf")

# Combine plots into a list
plot_list <- list(boxplot_tea, boxplot_wilcox)

# Save all plots with minimal page numbers
ggsave(output_pdf, marrangeGrob(plot_list, nrow=1, ncol=1, top=""),
       width = 12, height = 8, device = "pdf")













