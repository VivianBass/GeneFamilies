
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(dotenv)

output_data_dir <- Sys.getenv("OUTPUT_DATA_DIR")
results_dir <- Sys.getenv("RESULTS_DIR")


load(file.path(output_data_dir, "expression_profile_distances_statistics.RData"))

# Function to determine the level of significance
significance_level <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")  # Not significant
}


# combined dataframe

df_median_mean_paralogs, 
df_median_mean_orthologs,


# Create a dataframe with median distances
median_distances <- data.frame(
  distance = c(gene_family_expression_dists_stats$median_ortholog_exp_dists,
               gene_family_expression_dists_stats$median_paralog_exp_dists),
  group = rep(c("Orthologs", "Paralogs"), each = nrow(gene_family_expression_dists_stats))
)

# Create a dataframe with mean distances
mean_distances <- data.frame(
  distance = c(gene_family_expression_dists_stats$mean_ortholog_exp_dists,
               gene_family_expression_dists_stats$mean_paralog_exp_dists),
  group = rep(c("Orthologs", "Paralogs"), each = nrow(gene_family_expression_dists_stats))
)

# t-tests
# # Uncomment this part to perform t-test for medians
# t_test_medians <- t.test(distance ~ group, data = median_distances, alternative = 'greater')

# Uncomment this part to perform t-test for means
# t_test_mean_results <- t.test(distance ~ group, data = mean_distances, alternative = 'greater')

# Uncomment this part to save the t-test results in a dataframe
# t_test_summary <- data.frame(
#   analysis = c("Median", "Mean"),



  
#   # Results without correction
#   statistic_original = c(
#                           t_test_medians$statistic, 
#                           t_test_mean_results$statistic),
#   p_value_original = c(
#                        t_test_medians$p.value, 
#                        t_test_mean_results$p.value),
#   conf_int_lower_original = c(
#                                t_test_medians$conf.int[1], 
#                                t_test_mean_results$conf.int[1]),
#   conf_int_upper_original = c(
#                                t_test_medians$conf.int[2], 
#                                t_test_mean_results$conf.int[2]),
#   # Use the means reported by the t-tests
#   mean_orthologs = c(
#                       t_test_medians$estimate[1],     
#                       t_test_mean_results$estimate[1]),  
#   mean_paralogs = c(
#                      t_test_medians$estimate[2],     
#                      t_test_mean_results$estimate[2]), 
  
#   # Significance
#   significance_original = c(
#                     significance_level(t_test_medians$p.value),
#                     significance_level(t_test_mean_results$p.value)),
#   p_value_corrected = p.adjust(c(
#                                   t_test_medians$p.value, 
#                                   t_test_mean_results$p.value), 
#                                 method = "BH")
# )

# # Add a column of significance based on the corrected p-value
# t_test_summary$significance_corrected <- sapply(t_test_summary$p_value_corrected, significance_level)

# Create a boxplot for median distances
boxplot_median <- ggplot(median_distances, aes(x = group, y = distance, fill = group)) +
  geom_boxplot(outlier.shape = NA) +  # Exclude outliers from the boxplot to avoid overlap
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +  # Add jittered points
  labs(title = "Median Expression Distances",
       y = "Distance") +
  theme_minimal() +
  geom_signif(comparisons = list(c("Orthologs", "Paralogs")),
              map_signif_level = TRUE)

# Create a boxplot for mean distances
boxplot_mean <- ggplot(mean_distances, aes(x = group, y = distance, fill = group)) +
  geom_boxplot(outlier.shape = NA) +  # Exclude outliers from the boxplot to avoid overlap
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +  # Add jittered points
  labs(title = "Mean Expression Distances",
       y = "Distance") +
  theme_minimal() +
  geom_signif(comparisons = list(c("Orthologs", "Paralogs")),
              map_signif_level = TRUE)

# Save all three boxplots in a single PDF
pdf(file.path(results_dir,"boxplots_expression_distances_with_jitter.pdf"), height = 15, width = 7)
grid.arrange(boxplot_median, boxplot_mean, ncol = 1)
dev.off()

# # Uncomment this part to save the original and corrected t-test results to a CSV file
# write.csv(t_test_summary, file.path(results_dir, "t_test_results_original_and_corrected.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------



library(ggplot2)
library(tidyr)
library(dplyr)
library(ggthemes)



# Select relevant columns for plotting, 
col_names <- c("Mean_Fat_Body", "Mean_Gut", "Mean_Muscle", "Mean_Whole_Body")

#if you want to dynamically filter columns based on a pattern in their names (e.g., selecting all columns #that contain "count"), you could use the dplyr::select() with helper functions like contains():

col_names <- names(df %>% select(contains("count")))
col_names <- names(df %>% select_if(is.numeric))



plot_title <- "df_median/mean_exp.prof.dists.tissue"
x_axis_label <- ""
y_axis_label <- "exp.prof.dists_per_cluster"

# Pivot the dataframe to long format
df <- your_df

df_long <- df %>% pivot_longer(cols = all_of(col_names), names_to = "Category", values_to = "Count")

create_boxplot <- function(df_long, plot_title, y_axis_label, output_file = "median_mean_exp.prof.dists.tissue.pdf") {
  library(ggplot2)

  plot_vv <- ggplot(df_long, aes(x = Category, y = Count, fill = Category)) +
    geom_boxplot(color = "black", outlier.shape = 21, outlier.color = "red", outlier.fill = "red") +
    scale_fill_brewer(palette = "Pastel1") +
    labs(title = plot_title, x = " ", y = y_axis_label) +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.05)) +
    theme(
      plot.title = element_text(size = 14, face = "bold", margin = margin(t = 20, r = 0 ,b = 20, l = 0), hjust = 0.40),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 20, r = 0, b = 20, l = 0), hjust = 0.32),
      axis.title.y = element_text(size = 14, margin = margin(t = 20, r = 20, b = 20, l = 20)),
      axis.text.x = element_text(size = 10, angle = 35, hjust = 0.5, vjust = 0.5),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_line(color = "gray80"),
      legend.position = "none"
    )

  ggsave(output_file, plot = plot_vv, width = 8, height = 6)
  
  return(plot_vv)
}

-----------------------------------------------------------------------------------------

# ggplot 2 and ggrain, for Raincloud plots

library(ggplot2)
library(ggrain)
library(dplyr)
library(tidyr)


# Pivot the dataframe to long format

df <- plot_df_exp.prof.dists_ortho_para

df_long <- df %>% pivot_longer(cols = all_of(col_names), names_to = "Category", values_to = "Value")

# Select relevant columns for plotting
col_names <- c("paralogs", "orthologs")

# Define titles and labels for the plot
plot_title <- "Raincloud_Plot_of_Expression_Profile_Distances_para/ortho"
x_axis_label <- ""
y_axis_label <- "Expression Profile Distances"

create_raincloud_plot <- function(df_long, plot_title, x_axis_label, y_axis_label, output_file = "tissue.pdf") {

  library(ggplot2)
  library(ggrain)

  plot_vv <- ggplot(df_long, aes(x = Category, y = Value, fill = Category)) +
    geom_rain(alpha = 0.6, width = 0.6, justification = 0.5) +
    geom_boxplot(color = "black", width = 0.1, position = position_nudge(x = 0), alpha = 0.4) +
    geom_jitter(aes(color = Category), width = 0.2, size = 0.6, alpha = 0.4) +  
    theme_minimal(base_size = 14) +
    scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.05)) +
    scale_fill_brewer(palette = "Pastel1") +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = plot_title,
      x = x_axis_label,
      y = y_axis_label
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold", margin = margin(t = 10, b = 10), hjust = 0.5),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 15)),
      axis.title.y = element_text(size = 14, margin = margin(r = 10)),
      axis.text.x = element_text(size = 10, angle = 30, hjust = 0.5, vjust = 0.65),
      axis.text.y = element_text(size = 12),
      plot.margin = margin(t = 10, r = 30, b = 10, l = 10),
      panel.grid.major = element_line(color = "gray85"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )

  ggsave(output_file, plot = plot_vv, width = 8, height = 6)
  
  return(plot_vv)
}




