

# Function for Calculation of Means in a Dataframe

calculate_means <- function(df, column_name) {

  library(dplyr)

  df_mean <- df %>%
    rowwise() %>%  

    mutate(
      mean_value = mean(as.numeric(unlist(strsplit(as.character(!!sym(column_name)), ","))), na.rm = TRUE)
    ) %>%

    ungroup() %>%
    select(Family, mean_value)
  
  colnames(df_mean) <- c("Family", column_name)

  return(df_mean)
}

# -----------------------------------------------------------------------------

# Function for performing t-tests and Wilcoxon tests on pairwise distribution data

# Input: data (Dataframe) and output_file (Output file name)

perform_analysis <- function(data_long, output_file) {

  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(rstatix)
  library(ggpubr)  
  
  # Perform t-test and Wilcoxon test
  ttest_result <- t_test(value ~ group, data = data_long, paired = TRUE)
  wilcox_result <- wilcox_test(value ~ group, data = data_long, paired = TRUE)
  
  # Adjust p-values using Benjamini-Hochberg (BH) method
  p_vals <- c(ttest_result$p, wilcox_result$p)
  adjusted_p_vals <- p.adjust(p_vals, method = "BH")
  
  # Append adjusted p-values
  ttest_result <- ttest_result %>% mutate(adjusted_p = adjusted_p_vals[1:n()])
  wilcox_result <- wilcox_result %>% mutate(adjusted_p = adjusted_p_vals[(n() + 1):(2 * n())])
  
  # Combine results
  combined_results <- bind_rows(ttest_result, wilcox_result)
  
  # Define file path and write output to .tsv file
  write.table(combined_results, file = file.path("plots", paste0(output_file, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

  return(list(results = combined_results))
}

# Function for performing t-tests and Wilcoxon tests on pairwise distribution data

# Input: data (Dataframe) and output_file (Output file name)

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

# -----------------------------------------------------------------------------

# Function for creating and saving boxplot from t-test and Wilcoxon test results

# Input: data_long (Dataframe) and output_file (Output file name)

create_and_save_boxplot <- function(data_long, output_file = "boxplot_comparison.pdf",
                                    plot_title = NULL, x_label = NULL, y_label = NULL,
                                    custom_x_labels = NULL, point_size = 0.3,  # Reduced size for invisibility
                                    point_alpha = 0.3,  # Increase transparency
                                    plot_width = 14, plot_height = 8) {
  library(ggplot2)
  library(ggpubr)

  # Create the 'plots' directory if it doesn't exist
  if (!dir.exists("plots")) {
    dir.create("plots")
  }

  # Ensure 'group' is a factor for proper coloring
  data_long$group <- as.factor(data_long$group)

  # Create the boxplot with larger elements
  plot <- ggboxplot(data_long, x = "group", y = "value", fill = "group", 
                    add = "jitter", add.params = list(size = point_size, alpha = point_alpha)) +  # Adjust size and alpha
    scale_fill_brewer(palette = "Pastel1") +
    theme_minimal(base_size = 14) +  # Base size increased for all text elements
    scale_y_continuous(limits = c(min(data_long$value), max(data_long$value)), breaks = seq(0, 0.5, by = 0.05)) +
    theme(
      plot.title = element_text(size = 13, face = "bold", margin = margin(t = 20, r = 0 ,b = 20, l = 0), hjust = 0.40),
      axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 20, r = 0, b = 20, l = 0), hjust = 0.32),
      axis.title.y = element_text(size = 12, face = "bold", margin = margin(t = 20, r = 20, b = 20, l = 20)),
      axis.text.x = element_text(size = 10, angle = 35, hjust = 0.5, vjust = 0.5),
      axis.text.y = element_text(size = 10),
      plot.margin = margin(t = 10, r = 30, b = 10, l = 10),
      panel.grid.major = element_line(color = "gray80"),
      legend.position = "none"
    ) +
    # Increased label size and color for statistical comparisons
    stat_compare_means(method = "t.test", paired = TRUE, label = "p.signif", 
                      label.y = max(data_long$value) * 0.9, size = 10, color = "blue") +
    stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p.signif", 
                      label.y = max(data_long$value) * 0.75, size = 10, color = "red")

  # Add optional titles and labels
  if (!is.null(plot_title)) plot <- plot + ggtitle(plot_title)
  if (!is.null(x_label)) plot <- plot + xlab(x_label)
  if (!is.null(y_label)) plot <- plot + ylab(y_label)
  if (!is.null(custom_x_labels)) {
    plot <- plot + scale_x_discrete(labels = custom_x_labels)
  }

  # Save the plot to a larger PDF with increased dimensions
  ggsave(filename = file.path("plots", output_file), plot = plot, device = "pdf", width = plot_width, height = plot_height)

  return(plot)
}
