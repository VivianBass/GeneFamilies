

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

# Usage:
# df_means <- calculate_means(df = df, column_name = column_name)
# column_name <- "paralogs.exp.prof.dists"
# df <- para

# -----------------------------------------------------------------------------

# Function for performing t-tests and Wilcoxon tests on pairwise distribution data

# Input: data (Dataframe) and output_file (Output file name)

perform_analysis <- function(data, output_file) {

  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(rstatix)
  library(ggpubr)  
  
  ttest_result <- t_test(value ~ group, data = data_long, paired = TRUE)
  wilcox_result <- wilcox_test(value ~ group, data = data_long, paired = TRUE)
  
  p_vals <- c(ttest_result$p, wilcox_result$p)
  adjusted_p_vals <- p.adjust(p_vals, method = "BH")
  
  ttest_result <- ttest_result %>% mutate(adjusted_p = adjusted_p_vals[1:n()])
  wilcox_result <- wilcox_result %>% mutate(adjusted_p = adjusted_p_vals[(n() + 1):(2 * n())])
  
  combined_results <- bind_rows(ttest_result, wilcox_result)

  write.table(combined_results, file = paste0(output_file, ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  ggboxplot(data_long, x = "group", y = "value", add = "jitter") +
    stat_compare_means(method = "t.test", paired = TRUE, label = "p.signif") +
    stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p.signif")
  
  return(list(results = combined_results, plot = plot))
}

# Usage:
# result <- perform_analysis(your_dataframe, output_file = "results.txt")
# print(result$results)
# print(result$plot)

# -----------------------------------------------------------------------------

# Function for creating and saving boxplot from t-test and Wilcoxon test results

# Input: data_long (Dataframe) and output_file (Output file name)

create_and_save_boxplot <- function(data_long, output_file = "boxplot_comparison.pdf",
                                    plot_title = NULL, x_label = NULL, y_label = NULL) {
  library(ggplot2)
  library(ggpubr)

  plot <- ggboxplot(data_long, x = "group", y = "value", fill = "group", add = "jitter") +
    geom_boxplot(color = "black", outlier.shape = 21, outlier.color = "red", outlier.fill = "red") +
    scale_fill_brewer(palette = "Pastel1") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.05)) +
    theme(
      plot.title = element_text(size = 14, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.40),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 20), hjust = 0.32),
      axis.title.y = element_text(size = 14, margin = margin(r = 20)),
      axis.text.x = element_text(size = 10, angle = 35, hjust = 0.5),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_line(color = "gray80"),
      legend.position = "none"
    ) +
    stat_compare_means(method = "t.test", paired = TRUE, label = "p.format") +
    stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p.format", label.y = 0.45)

  if (!is.null(plot_title)) plot <- plot + ggtitle(plot_title)
  if (!is.null(x_label)) plot <- plot + xlab(x_label)
  if (!is.null(y_label)) plot <- plot + ylab(y_label)

  ggsave(filename = output_file, plot = plot, device = "pdf", width = 8, height = 6)
 
  return(plot)
}

# Usage:

#data <- para_ortho_mean_df_mean
#cols <- colnames(data)
#cols <- cols[cols != "Family"]

#plot <- create_and_save_boxplot(
#  data_long = data_long,  
#  output_file = output_file,
#  plot_title = plot_title,
#  x_label = x_label,
#  y_label = y_label
#)