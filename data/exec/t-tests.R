# 1. T-Test (t.test): Used to determine if the means of two groups (e.g., orthologs vs paralogs) are significantly different. For example, you can compare mean median expression profile distances between groups using t.test(). The alternative hypothesis (e.g., "greater", "less") specifies the direction of the test.

# 2. Wilcox on Rank Sum Test (wilcox.test): A non-parametric test that checks if the overall distributions between two groups differ without assuming normality. It assesses whether one group tends to have larger values than the other.

# 3. P-Value Adjustment (p.adjust): After obtaining p-values from multiple comparisons (t-tests and Wilcoxon tests), they are adjusted to control for multiple hypothesis testing using p.adjust(). This helps to reduce false positives.

# Benjamini-Hochberg (BH) Correction: This method is used to adjust p-values for multiple hypothesis testing by controlling the False Discovery Rate (FDR). It ensures a balance between detecting true effects and minimizing false discoveries.

# -----------------------------------------------------------------------------

# R Libraries for t-tests and Wilcoxon tests, rstatix and ggpubr

# -----------------------------------------------------------------------------

# t-tests and wilcox-test on pairwise distribution data

1. Distributions of distances of paralogs and orthologs
  2. Distribution of mean distances
  3. Distribution of median distances

4. Distributions of distances per tissue
  5. Distribution of mean distances per tissue
  6. Distribution of median distances per tissue

7. Expression angles
8. Expression versatility

- example t_test_examples.R

a <- rnorm(10000, 100, 5)
b <- rnorm(10000, 1000, 50)
summary(a)
summary(b)
boxplot(list(a=a, b=b))
t.test(a, b, alternative='greater')
t.test(b, a, alternative='greater')
?t.test
t.test(b, a, alternative='greater')
savehistory('~/Desktop/t_test_examples.R')




load("y-data.RData")

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
df <- para_mean

library(dplyr)

df_means <- calculate_means(df , column_name = column_name)

para_mean_df_mean <- df_means
para_mean_df_median <- df_means
ortho_mean_df_mean <- df_means
ortho_mean_df_median <- df_means

# merging the dataframes on which the pairwise analysis is performed

merged_df <- ortho_mean_df_median %>%
  left_join(para_mean_df_median, by = "Family")

para_ortho_mean_df_mean <- merged_df
para_ortho_mean_df_median <- merged_df

# performing t-tests and Wilcoxon tests on pairwise distribution data

data <- para_ortho_mean_df_mean

result <- perform_analysis(data = data, output_file = "results")


# provide the output file name, Achsis labels etc. as arguments

data <- para_ortho_mean_df_mean
cols <- colnames(data)
cols <- cols[cols != "Family"]

data_long <- data %>% pivot_longer(cols = all_of(cols), names_to = "group", values_to = "value")

plot_title <- "Boxplot Comparison"
x_label <- " "
y_label <- "Mean Expression Profile Distance"
output_file <- "boxplot_comparison_result.pdf"

plot <- create_and_save_boxplot(
  data_long = data_long,  
  output_file = output_file,
  plot_title = plot_title,
  x_label = x_label,
  y_label = y_label
)

# ------------------------------------------------------------------------------------------------------
4. Distributions of distances per tissue
# ------------------------------------------------------------------------------------------------------

  5. Distribution of mean distances per tissue
  6. Distribution of median distances per tissue

  df_median_mean_paralogs.exp.prof.dists.tissue
  df_median_mean_orthologs.exp.prof.dists.tissue

# adjusting and selecting the columns of interest for the pairwise analysis

df_ortho_mean <- df_median_mean_paralogs.exp.prof.dists.tissue[ ,  c("Cluster", "Mean_Fat_Body", "Mean_Gut", "Mean_Muscle", "Mean_Whole_Body")]
df_ortho_median <- df_median_mean_orthologs.exp.prof.dists.tissue[ ,  c("Cluster", "Median_Fat_Body", "Median_Gut", "Median_Muscle", "Median_Whole_Body")]
df_para_mean <- df_median_mean_paralogs.exp.prof.dists.tissue[ ,  c("Cluster", "Mean_Fat_Body", "Mean_Gut", "Mean_Muscle", "Mean_Whole_Body")]
df_para_median <- df_median_mean_orthologs.exp.prof.dists.tissue[ ,  c("Cluster", "Median_Fat_Body", "Median_Gut", "Median_Muscle", "Median_Whole_Body")]

# merging the dataframes on which the pairwise analysis is performed

merged_df <- df_ortho_median %>%
  left_join(df_para_median, by = "Cluster")

df_ortho_para_mean <- merged_df
df_ortho_para_median <- merged_df

# conerting the dataframes to long format for t-tests and Wilcoxon tests

data <- df_ortho_para_median
cols <- colnames(data)
cols <- cols[cols != "Cluster"]

data_long <- data %>%
  pivot_longer(cols = all_of(cols), 
               names_to = "group", 
               values_to = "value")

source("data/exec/t_test_functions.R")

# performing t-tests and Wilcoxon tests on pairwise distribution data
# stored in the tsv file

result <- perform_analysis(data = data, output_file = "results-df_ortho_para_median")

# Plotting the t-test results

plot_title <- "Boxplot Comparison"
x_label <- " "
y_label <- "Mean Expression Profile Distance"
output_file <- "boxplot_comparison_result_jdfnjdfn.pdf"

plot <- create_and_save_boxplot(
  data_long = data_long,  
  output_file = output_file,
  plot_title = plot_title,
  x_label = x_label,
  y_label = y_label
)
