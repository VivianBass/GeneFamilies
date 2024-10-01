# 1. T-Test (t.test): Used to determine if the means of two groups (e.g., orthologs vs paralogs) are significantly different. For example, you can compare mean median expression profile distances between groups using t.test(). The alternative hypothesis (e.g., "greater", "less") specifies the direction of the test.

# 2. Wilcox on Rank Sum Test (wilcox.test): A non-parametric test that checks if the overall distributions between two groups differ without assuming normality. It assesses whether one group tends to have larger values than the other.

# 3. P-Value Adjustment (p.adjust): After obtaining p-values from multiple comparisons (t-tests and Wilcoxon tests), they are adjusted to control for multiple hypothesis testing using p.adjust(). This helps to reduce false positives.

# Benjamini-Hochberg (BH) Correction: This method is used to adjust p-values for multiple hypothesis testing by controlling the False Discovery Rate (FDR). It ensures a balance between detecting true effects and minimizing false discoveries.

# -----------------------------------------------------------------------------

# R Libraries for Wilcoxon tests, rstatix and ggpubr

perform_wilcox_test <- function(data_long, output_file) {

  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(rstatix)
  library(ggpubr)  
  
  # Perform the Wilcoxon test
  wilcox_result <- data_long %>%
    wilcox_test(value ~ group, paired = TRUE) %>%  # If your data is paired
    adjust_pvalue(method = "BH") %>%               # Adjust p-values using Benjamini-Hochberg (BH)
    add_significance()                             # Add significance stars
  
  # Save the results to a .tsv file
  write.table(wilcox_result, 
              file = file.path("plots", paste0(output_file, ".tsv")), 
              sep = "\t", row.names = FALSE, quote = FALSE)

  # Return results for further analysis or plotting
  return(list(results = wilcox_result))
}

# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(rstatix)
library(ggpubr)

# Reshape the data to long format
mydata <- p.df  
mydata_long <- mydata %>% 
  pivot_longer(cols = -gene.type, names_to = "variables", values_to = "value") %>%
  filter(!is.na(value) & !is.infinite(value))

# Perform Wilcoxon test with Benjamini-Hochberg (BH) p-value adjustment
stat.test <- mydata_long %>%
  group_by(variables) %>%
  wilcox_test(value ~ gene.type, p.adjust.method = "BH") %>%
  add_significance()

# Create boxplot for the data
myplot <- ggboxplot(
  mydata_long,
  x = "gene.type",
  y = "value",
  fill = "gene.type",
  palette = "npg",
  legend = "none",
  ggtheme = theme_pubr(border = TRUE)
) +
  facet_wrap(~variables)

# Add p-values to the plot
stat.test <- stat.test %>% add_xy_position(x = "gene.type")
myplot_with_pvals <- myplot + stat_pvalue_manual(stat.test, label = "p.adj.signif")

# Save the plot
ggsave("plots/wilcox-test-plots/order_tissues_mean_boxplot_with_pvalues.pdf", plot = myplot_with_pvals, width = 12, height = 8)

# Save the statistical test results to a file
write_tsv(stat.test, "plots/wilcox-test-plots/order_tissues_mean_statistical_results.tsv")

# -----------------------------------------------------------------------------
1. Distributions of distances of paralogs and orthologs
# -----------------------------------------------------------------------------

  2. Distribution of mean distances
  3. Distribution of median distances

  - Data used:
    orthologs.exp.prof.dists_stats_df
    paralogs.exp.prof.dists_stats_df

# adjusting and selecting the columns of interest for the pairwise analysis

para_mean <- paralogs.exp.prof.dists_stats_df[ ,  c("Family", "mean")]
para_median <- paralogs.exp.prof.dists_stats_df[ ,  c("Family", "median")]

ortho_mean <- orthologs.exp.prof.dists_stats_df[ ,  c("Family", "mean")]
ortho_median <- orthologs.exp.prof.dists_stats_df[ ,  c("Family", "median")]

# calculate mean Values for the Dataframe if row/column Element contain multiple values
# column_name , column where the mean is calculated

column_name <- "mean"
df <- ortho_mean

library(dplyr)

source("R/t_test_functions.R")

df_means <- calculate_means(df , column_name = column_name)

para_mean_df_mean <- df_means
para_mean_df_median <- df_means
ortho_mean_df_mean <- df_means
ortho_mean_df_median <- df_means

# merging the dataframes on which the pairwise analysis is performed
# performing t-tests and Wilcoxon tests on pairwise distribution data
# provide the output file name, Achsis labels etc. as arguments to the function

merged_df <- ortho_mean_df_mean %>%
  left_join(para_mean_df_mean, by = "Family")
para_ortho_mean_df_mean <- merged_df

data <- para_ortho_mean_df_mean
output_file_2 <- "2._Distribution_of_mean_distances_(para_ortho)_(t-test_wilcox-test)"
cols <- colnames(data)
cols <- cols[cols != "Family"]
data_long <- data %>% pivot_longer(cols = all_of(cols), names_to = "group", values_to = "value")

source("exec/t_test_functions.R")
result <- perform_analysis(data_long, output_file_2)

plot_title <- output_file_2
x_label <- " "
y_label <- "mean_Expression_Profile_Distances"
custom_x_labels = c(
  "mean.x" = "Orthologs",
  "mean.y" = "Paralogs"
)
output_file <- file.path("plots", paste0(output_file_2, ".pdf"))

plot <- create_and_save_boxplot(
  data_long = data_long,  
  output_file = output_file,
  plot_title = plot_title,
  x_label = x_label,
  y_label = y_label,
  custom_x_labels = custom_x_labels,
  point_size = 0.5, 
  plot_width = 8, 
  plot_height = 6 
)


# -----------------------------------------------------------------------------

merged_df <- ortho_mean_df_median %>%
  left_join(para_mean_df_median, by = "Family")
para_ortho_mean_df_median <- merged_df

data <- para_ortho_mean_df_median
output_file_3 <- "3._Distribution_of_median_distances_(para_ortho)_(t-test_wilcox-test)"
cols <- colnames(data)
cols <- cols[cols != "Family"]
data_long <- data %>% pivot_longer(cols = all_of(cols), names_to = "group", values_to = "value")

source("exec/t_test_functions.R")
result <- perform_analysis(data_long, output_file_3)

plot_title <- output_file_3
x_label <- " "
y_label <- "median_Expression_Profile_Distances"
custom_x_labels = c(
  "median.x" = "Orthologs",
  "median.y" = "Paralogs"
)
output_file <- file.path("plots", paste0(output_file_3, ".pdf"))

plot <- create_and_save_boxplot(
  data_long = data_long,  
  output_file = output_file,
  plot_title = plot_title,
  x_label = x_label,
  y_label = y_label,
  custom_x_labels = custom_x_labels,
  point_size = 0.5, 
  plot_width = 8, 
  plot_height = 6 
)



# ------------------------------------------------------------------------------------------------------
4. Distributions of distances per tissue
# ------------------------------------------------------------------------------------------------------

  5. Distribution of mean distances per tissue
  6. Distribution of median distances per tissue

  df_median_mean_paralogs.exp.prof.dists.tissue
  df_median_mean_orthologs.exp.prof.dists.tissue

# adjusting and selecting the columns of interest for the pairwise analysis

df_ortho_mean_tissue <- df_median_mean_paralogs.exp.prof.dists.tissue[ ,  c("Cluster", "Mean_Fat_Body", "Mean_Gut", "Mean_Muscle", "Mean_Whole_Body")]
df_ortho_median_tissue <- df_median_mean_orthologs.exp.prof.dists.tissue[ ,  c("Cluster", "Median_Fat_Body", "Median_Gut", "Median_Muscle", "Median_Whole_Body")]
df_para_mean_tissue <- df_median_mean_paralogs.exp.prof.dists.tissue[ ,  c("Cluster", "Mean_Fat_Body", "Mean_Gut", "Mean_Muscle", "Mean_Whole_Body")]
df_para_median_tissue <- df_median_mean_orthologs.exp.prof.dists.tissue[ ,  c("Cluster", "Median_Fat_Body", "Median_Gut", "Median_Muscle", "Median_Whole_Body")]


# merging the dataframes on which the pairwise analysis is performed
# converting the dataframes to long format for t-tests and Wilcoxon tests, and performing Analysis
# Plotting the t-test results

source("data/exec/t_test_functions.R")

merged_df <- df_ortho_mean_tissue %>%
  left_join(df_para_mean_tissue, by = "Cluster")
df_ortho_para_mean_tissue <- merged_df

data <- df_ortho_para_mean_tissue
cols <- colnames(data)
cols <- cols[cols != "Cluster"]
data_long <- data %>% pivot_longer(cols = all_of(cols), names_to = "group", values_to = "value")
output_file_5 <- "5._Distribution_of_mean_distances_per_tissue_(para_ortho)_(t-test_wilcox-test)"
result <- perform_analysis(data_long, output_file_5)

plot_title <- output_file_5
x_label <- " "
y_label <- "mean_Expression_Profile_Distances"
custom_x_labels = c(
  "Mean_Fat_Body.x" = "Fat Body (Orthologs)",
  "Mean_Gut.x" = "Gut (Orthologs)",
  "Mean_Muscle.x" = "Muscle (Orthologs)",
  "Mean_Whole_Body.x" = "Whole Body (Orthologs)",
  "Mean_Fat_Body.y" = "Fat Body (Paralogs)",
  "Mean_Gut.y" = "Gut (Paralogs)",
  "Mean_Muscle.y" = "Muscle (Paralogs)",
  "Mean_Whole_Body.y" = "Whole Body (Paralogs)"
)

output_file <- file.path("plots", paste0(output_file_5, ".pdf"))

plot <- create_and_save_boxplot(
  data_long = data_long,  
  output_file = output_file,
  plot_title = plot_title,
  x_label = x_label,
  y_label = y_label,
  custom_x_labels = custom_x_labels,
  point_size = 0.1, 
  plot_width = 14, 
  plot_height = 8 
)


# -------------------------------------------------------

merged_df <- df_ortho_median_tissue %>%
  left_join(df_para_median_tissue, by = "Cluster")
df_ortho_para_median_tissue <- merged_df

data <- df_ortho_para_median_tissue
cols <- colnames(data)
cols <- cols[cols != "Cluster"]
data_long <- data %>% pivot_longer(cols = all_of(cols), names_to = "group", values_to = "value")
output_file_6 <- "6._Distribution_of_median_distances_per_tissue_(para_ortho)_(t-test_wilcox-test)"
result <- perform_analysis(data_long, output_file_6)

plot_title <- output_file_6
x_label <- " "
y_label <- "median_Expression_Profile_Distances"
custom_x_labels = c(
  "Median_Fat_Body.x" = "Fat Body (Orthologs)",
  "Median_Gut.x" = "Gut (Orthologs)",
  "Median_Muscle.x" = "Muscle (Orthologs)",
  "Median_Whole_Body.x" = "Whole Body (Orthologs)",
  "Median_Fat_Body.y" = "Fat Body (Paralogs)",
  "Median_Gut.y" = "Gut (Paralogs)",
  "Median_Muscle.y" = "Muscle (Paralogs)",
  "Median_Whole_Body.y" = "Whole Body (Paralogs)"
)

output_file <- file.path("plots", paste0(output_file_6, ".pdf"))

plot <- create_and_save_boxplot(
  data_long = data_long,  
  output_file = output_file,
  plot_title = plot_title,
  x_label = x_label,
  y_label = y_label,
  custom_x_labels = custom_x_labels,
  point_size = 0.1, 
  plot_width = 14, 
  plot_height = 8 
)

# ------------------------------------------------------------------------------------------------------
7. Distributions of Angles
# ------------------------------------------------------------------------------------------------------

8. Expression angles
9 Expression versatility