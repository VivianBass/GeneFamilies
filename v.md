
## remove redundant Code



```R
create_distance_boxplot <- function(data,
                                  type_combinations,
                                  title_prefix = "Expression Distances",
                                  test_type = "t.test",
                                  alternative = "two.sided",
                                  plot_width = 12,
                                  plot_height = 8,
                                  nrow = 2,
                                  ncol = 2) {
    
    ggplot(data, aes(x = Type, y = Distance, fill = Type)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
        labs(title = paste(title_prefix, paste0("(", test_type, ")")), y = "Distance") +
        theme_pubr(border = TRUE) +
        scale_y_continuous(breaks = seq(0, max(data$Distance), by = 0.2)) +
        geom_signif(
            comparisons = type_combinations,
            test = test_type,
            test.args = list(alternative = alternative),
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
}
```

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------


```R
create_and_save_distance_boxplots <- function(mean_data,
                                            median_data,
                                            type_combinations,
                                            alternative = "two.sided",
                                            results_dir,
                                            filename = "boxplots_expression_distances_all_tests.pdf",
                                            plot_width = 16,
                                            plot_height = 12,
                                            nrow = 2,
                                            ncol = 2) {
    
    # T-test plots
    boxplot_mean_ttest <- create_distance_boxplot(
        data = mean_data,
        type_combinations = type_combinations,
        title_prefix = "Mean Expression Distances",
        test_type = "t.test",
        alternative = alternative,
        plot_width = plot_width,
        plot_height = plot_height,
        nrow = nrow,
        ncol = ncol
    )
    
    boxplot_median_ttest <- create_distance_boxplot(
        data = median_data,
        type_combinations = type_combinations,
        title_prefix = "Median Expression Distances",
        test_type = "t.test",
        alternative = alternative,
        plot_width = plot_width,
        plot_height = plot_height,
        nrow = nrow,
        ncol = ncol
    )
    
    # Wilcox test plots
    boxplot_mean_wilcox <- create_distance_boxplot(
        data = mean_data,
        type_combinations = type_combinations,
        title_prefix = "Mean Expression Distances",
        test_type = "wilcox.test",
        alternative = alternative,
        plot_width = plot_width,
        plot_height = plot_height,
        nrow = nrow,
        ncol = ncol
    )
    
    boxplot_median_wilcox <- create_distance_boxplot(
        data = median_data,
        type_combinations = type_combinations,
        title_prefix = "Median Expression Distances",
        test_type = "wilcox.test",
        alternative = alternative,
        plot_width = plot_width,
        plot_height = plot_height,
        nrow = nrow,
        ncol = ncol
    )
    
    # Combine all plots and save
    output_pdf <- file.path(results_dir, filename)
    plot_list <- list(boxplot_median_ttest, boxplot_mean_ttest,
                     boxplot_median_wilcox, boxplot_mean_wilcox)
    ggsave(output_pdf,
           marrangeGrob(plot_list, nrow=nrow, ncol=ncol, top=""),
           width = plot_width, height = plot_height,
           device = "pdf")
    
    return(plot_list)
}
```

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------


```R
process_distance_statistics <- function(input_file, stats_pattern, output_file, output_data_dir) {
    # Load data
    load(file.path(output_data_dir, input_file))
    loaded_objects <- ls()
    data_names <- loaded_objects[grepl(stats_pattern, loaded_objects)]
    
    # Initialize empty data frames
    df_mean.dists <- data.frame()
    df_median.dists <- data.frame()
    
    # Process datasets
    for (data_name in data_names) {
        current_data <- get(data_name)
        type_name <- sub(stats_pattern, "", data_name)
        
        # Compile mean distances
        temp_mean_df <- tibble(
            Type = type_name,
            Cluster = names(current_data$Mean),
            Distance = unlist(current_data$Mean)
        )
        df_mean.dists <- bind_rows(df_mean.dists, temp_mean_df)
        
        # Compile median distances
        temp_median_df <- tibble(
            Type = type_name,
            Cluster = names(current_data$Median),
            Distance = unlist(current_data$Median)
        )
        df_median.dists <- bind_rows(df_median.dists, temp_median_df)
    }
    
    # Save processed dataframes
    save(df_mean.dists, df_median.dists,
         file = file.path(output_data_dir, output_file))
    
    return(list(mean = df_mean.dists, median = df_median.dists))
}
```


# ---------------------------------------------------------------------
# ---------------------------------------------------------------------


```R
# Process regular distance statistics
regular_dfs <- process_distance_statistics(
    input_file = "exp.prof.dists_statistics.RData",
    stats_pattern = ".lst_dists_stats$",
    output_file = "exp.prof.dists_mean_median.RData",
    output_data_dir = output_data_dir
)

# Process log2 distance statistics
log2_dfs <- process_distance_statistics(
    input_file = "exp.prof.dists_statistics_log2.RData",
    stats_pattern = ".lst_dists_log2_stats$",
    output_file = "exp.prof.dists_mean_median_log2.RData",
    output_data_dir = output_data_dir
)

# Process angles distance statistics
angles_dfs <- process_distance_statistics(
    input_file = "exp.prof.angels_statistics.RData",
    stats_pattern = ".lst_cos_angles_dists_stats$",
    output_file = "exp.prof.angles_mean_median.RData",
    output_data_dir = output_data_dir
)

# Process angles log2 distance statistics
angles_log2_dfs <- process_distance_statistics(
    input_file = "exp.prof.angels_statistics_log2.RData",
    stats_pattern = ".lst_cos_angles_dists_log2_stats$",
    output_file = "exp.prof.angles_mean_median_log2.RData",
    output_data_dir = output_data_dir
)

# ---------------------------------------------------------------------

# Create plots for regular data
types <- unique(regular_dfs$mean$Type)
type_combinations <- combn(types, 2, simplify = FALSE)
plots_regular <- create_and_save_distance_boxplots(
    mean_data = regular_dfs$mean,
    median_data = regular_dfs$median,
    type_combinations = type_combinations,
    results_dir = results_dir,
    filename = "boxplots_expression_distances_regular_all_tests.pdf"
)

# Create plots for log2 data
types <- unique(log2_dfs$mean$Type)
type_combinations <- combn(types, 2, simplify = FALSE)
plots_log2 <- create_and_save_distance_boxplots(
    mean_data = log2_dfs$mean,
    median_data = log2_dfs$median,
    type_combinations = type_combinations,
    results_dir = results_dir,
    filename = "boxplots_expression_distances_log2_all_tests.pdf"
)

# Create plots for angles data
types <- unique(angles_dfs$mean$Type)
type_combinations <- combn(types, 2, simplify = FALSE)
plots_angles <- create_and_save_distance_boxplots(
    mean_data = angles_dfs$mean,
    median_data = angles_dfs$median,
    type_combinations = type_combinations,
    results_dir = results_dir,
    filename = "boxplots_expression_distances_angles_all_tests.pdf"
)

# Create plots for angles log2 data
types <- unique(angles_log2_dfs$mean$Type)
type_combinations <- combn(types, 2, simplify = FALSE)
plots_angles_log2 <- create_and_save_distance_boxplots(
    mean_data = angles_log2_dfs$mean,
    median_data = angles_log2_dfs$median,
    type_combinations = type_combinations,
    results_dir = results_dir,
    filename = "boxplots_expression_distances_angles_log2_all_tests.pdf"
)
```