
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

# Filter out rows with NaN values in angle.diag
plot.df <- p.df[!is.nan(p.df$angle.diag), ]

# Set the factor levels for gene types
plot.df$gene.type <- factor(plot.df$gene.type, levels = unique(plot.df$gene.type))

# Define color palette
colors <- brewer.pal(length(unique(plot.df$gene.type)), "Pastel1")

# --------------------------------------------------------------------------

# Create gene type combinations first
gene_types <- levels(plot.df$gene.type)
type_combinations <- combn(gene_types, 2, simplify = FALSE)

# Create the angle plot with enhanced significance annotations
ggplot_angle <- ggplot(plot.df, aes(x = gene.type, y = angle.diag, fill = gene.type)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "Expression Angle To Diagonal",
       y = "Relative Tissue Specificity",
       x = "Gene Type") +
  theme_pubr(border = TRUE) +
  scale_y_continuous(breaks = seq(0, max(plot.df$angle.diag, na.rm = TRUE), by = 0.1)) +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
    axis.text.x = element_text(size = 10),
    plot.margin = margin(r = 30)
  ) +
  geom_signif(
    comparisons = type_combinations,
    test = "t.test",
    test.args = list(alternative = "greater"),
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
    step_increase = 0.05,
    tip_length = 0.005,
    vjust = 0.5,
    color = "black",
    size = 0.3,
    textsize = 2.5
  )

# Save the angle plot
ggsave(file.path(results_dir, "expressionAngleToDiagonalBoxplot_t.test.pdf"),
       ggplot_angle, width = 10, height = 8)


# Create the angle plot with enhanced significance annotations
ggplot_angle <- ggplot(plot.df, aes(x = gene.type, y = angle.diag, fill = gene.type)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "Expression Angle To Diagonal",
       y = "Relative Tissue Specificity",
       x = "Gene Type") +
  theme_pubr(border = TRUE) +
  scale_y_continuous(breaks = seq(0, max(plot.df$angle.diag, na.rm = TRUE), by = 0.1)) +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
    axis.text.x = element_text(size = 10),
    plot.margin = margin(r = 30)
  ) +
  geom_signif(
    comparisons = type_combinations,
    test = "wilcox.test",
    test.args = list(alternative = "greater"),
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
    step_increase = 0.05,
    tip_length = 0.005,
    vjust = 0.5,
    color = "black",
    size = 0.3,
    textsize = 2.5
  )

# Save the angle plot
ggsave(file.path(results_dir, "expressionAngleToDiagonalBoxplot_wilcox.test.pdf"),
       ggplot_angle, width = 10, height = 8)

# --------------------------------------------------------------------------

# Calculate relative versatility and add as a new column
plot.df$rel.vers <- 1 - plot.df$angle.diag

# Create the versatility plot with enhanced significance annotations
ggplot_vers <- ggplot(plot.df, aes(x = gene.type, y = rel.vers, fill = gene.type)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "Relative Expression Versatility",
       y = "Relative Tissue Versatility",
       x = "Gene Type") +
  theme_pubr(border = TRUE) +
  scale_y_continuous(breaks = seq(0, max(plot.df$rel.vers, na.rm = TRUE), by = 0.1)) +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
    axis.text.x = element_text(size = 10),
    plot.margin = margin(r = 30)
  ) +
  geom_signif(
    comparisons = type_combinations,
    test = "t.test",
    test.args = list(alternative = "two.sided"),
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
    step_increase = 0.05,
    tip_length = 0.005,
    vjust = 0.5,
    color = "black",
    size = 0.3,
    textsize = 2.5
  )

# Save the versatility plot
ggsave(file.path(results_dir, "relativeExpressionVersatilityBoxplot_t.test.pdf"),
       ggplot_vers, width = 10, height = 8)


ggplot_vers <- ggplot(plot.df, aes(x = gene.type, y = rel.vers, fill = gene.type)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "Relative Expression Versatility",
       y = "Relative Tissue Versatility",
       x = "Gene Type") +
  theme_pubr(border = TRUE) +
  scale_y_continuous(breaks = seq(0, max(plot.df$rel.vers, na.rm = TRUE), by = 0.1)) +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
    axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
    axis.text.x = element_text(size = 10),
    plot.margin = margin(r = 30)
  ) +
  geom_signif(
    comparisons = type_combinations,
    test = "wilcox.test",
    test.args = list(alternative = "greater"),
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
    step_increase = 0.05,
    tip_length = 0.005,
    vjust = 0.5,
    color = "black",
    size = 0.3,
    textsize = 2.5
  )

# Save the versatility plot
ggsave(file.path(results_dir, "relativeExpressionVersatilityBoxplot_wilcox.test.pdf"),
       ggplot_vers, width = 10, height = 8)




