





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






process_regular_angles <- function(data_pattern, loaded_objects) {
    data_names <- loaded_objects[grepl(data_pattern, loaded_objects)]
    
    data_names %>%
        map_df(function(data_name) {
            current_data <- get(data_name)
            type_name <- sub(data_pattern, "", data_name)
            
            tibble(
                Type = type_name,
                Angle = unlist(current_data)
            )
        })
}


create_angle_boxplot <- function(data, type_combinations, test_type, y_breaks, title_prefix = "Complete Expression Angles") {
    ggplot(data, aes(x = Type, y = Angle, fill = Type)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
        labs(title = paste(title_prefix, paste0("(", test_type, ")")), y = "Angle") +
        theme_pubr(border = TRUE) +
        scale_y_continuous(breaks = y_breaks) +
        geom_signif(
            comparisons = type_combinations,
            test = test_type,
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
}

plot_and_save_combined_angles <- function(regular_data, log2_data, results_dir, filename_prefix) {
    # Process regular data
    df_regular_filtered <- regular_data %>% filter(is.finite(Angle))
    types_regular <- unique(df_regular_filtered$Type)
    type_combinations_regular <- combn(types_regular, 2, simplify = FALSE)
    y_breaks_regular <- seq(0, max(df_regular_filtered$Angle, na.rm = TRUE), by = 0.2)
    
    # Process log2 data
    df_log2_filtered <- log2_data %>% filter(is.finite(Angle))
    types_log2 <- unique(df_log2_filtered$Type)
    type_combinations_log2 <- combn(types_log2, 2, simplify = FALSE)
    y_breaks_log2 <- seq(0, max(df_log2_filtered$Angle, na.rm = TRUE), by = 0.2)
    
    # Create all plots
    plot_list <- list(
        create_angle_boxplot(df_regular_filtered, type_combinations_regular, "t.test", y_breaks_regular, "Regular Expression Angles"),
        create_angle_boxplot(df_log2_filtered, type_combinations_log2, "t.test", y_breaks_log2, "Log2 Expression Angles"),
        create_angle_boxplot(df_regular_filtered, type_combinations_regular, "wilcox.test", y_breaks_regular, "Regular Expression Angles"),
        create_angle_boxplot(df_log2_filtered, type_combinations_log2, "wilcox.test", y_breaks_log2, "Log2 Expression Angles")
    )
    
    # Save combined plots
    output_pdf <- file.path(results_dir, paste0(filename_prefix, "_combined.pdf"))
    ggsave(output_pdf,
           marrangeGrob(plot_list, nrow=2, ncol=2, top=""),
           width = 16, height = 12,
           device = "pdf")
    
    return(plot_list)
}





create_boxplot <- function(data, type_combinations, test_type, y_breaks, title_prefix = "Complete Expression Distances") {
    ggplot(data, aes(x = Type, y = Distance, fill = Type)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
        labs(title = paste(title_prefix, paste0("(", test_type, ")")), y = "Distance") +
        theme_pubr(border = TRUE) +
        scale_y_continuous(breaks = y_breaks) +
        geom_signif(
            comparisons = type_combinations,
            test = test_type,
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
}




plot_and_save_combined_distances <- function(regular_data, log2_data, results_dir, filename_prefix) {
    # Process regular data
    df_regular_filtered <- regular_data %>% filter(is.finite(Distance))
    types_regular <- unique(df_regular_filtered$Type)
    type_combinations_regular <- combn(types_regular, 2, simplify = FALSE)
    y_breaks_regular <- seq(0, max(df_regular_filtered$Distance, na.rm = TRUE), by = 0.2)
    
    # Process log2 data
    df_log2_filtered <- log2_data %>% filter(is.finite(Distance))
    types_log2 <- unique(df_log2_filtered$Type)
    type_combinations_log2 <- combn(types_log2, 2, simplify = FALSE)
    y_breaks_log2 <- seq(0, max(df_log2_filtered$Distance, na.rm = TRUE), by = 0.2)
    
    # Create all plots
    plot_list <- list(
        create_boxplot(df_regular_filtered, type_combinations_regular, "t.test", y_breaks_regular, "Regular Expression Distances"),
        create_boxplot(df_log2_filtered, type_combinations_log2, "t.test", y_breaks_log2, "Log2 Expression Distances"),
        create_boxplot(df_regular_filtered, type_combinations_regular, "wilcox.test", y_breaks_regular, "Regular Expression Distances"),
        create_boxplot(df_log2_filtered, type_combinations_log2, "wilcox.test", y_breaks_log2, "Log2 Expression Distances")
    )
    
    # Save combined plots
    output_pdf <- file.path(results_dir, paste0(filename_prefix, "_combined.pdf"))
    ggsave(output_pdf,
           marrangeGrob(plot_list, nrow=2, ncol=2, top=""),
           width = 16, height = 12,
           device = "pdf")
    
    return(plot_list)
}




plot_and_save_combined_tissue_distances <- function(regular_data, log2_data, results_dir, filename_prefix) {
    # Filter data
    df_regular_filtered <- regular_data %>% filter(is.finite(Distance))
    df_log2_filtered <- log2_data %>% filter(is.finite(Distance))
    
    # Get unique tissues
    tissue_types <- unique(df_regular_filtered$Tissue)
    
    # Create plots for each tissue
    for (tissue in tissue_types) {
        df_regular_tissue <- subset(df_regular_filtered, Tissue == tissue)
        df_log2_tissue <- subset(df_log2_filtered, Tissue == tissue)
        
        types_regular <- unique(df_regular_tissue$Type)
        types_log2 <- unique(df_log2_tissue$Type)
        
        if (length(types_regular) >= 2 && length(types_log2) >= 2) {
            type_combinations_regular <- combn(types_regular, 2, simplify = FALSE)
            type_combinations_log2 <- combn(types_log2, 2, simplify = FALSE)
            
            y_breaks_regular <- seq(0, max(df_regular_tissue$Distance, na.rm = TRUE), by = 0.2)
            y_breaks_log2 <- seq(0, max(df_log2_tissue$Distance, na.rm = TRUE), by = 0.2)
            
            plot_list <- list(
                create_boxplot(df_regular_tissue, type_combinations_regular, "t.test", y_breaks_regular, 
                             paste("Regular Expression Distances -", tissue)),
                create_boxplot(df_log2_tissue, type_combinations_log2, "t.test", y_breaks_log2, 
                             paste("Log2 Expression Distances -", tissue)),
                create_boxplot(df_regular_tissue, type_combinations_regular, "wilcox.test", y_breaks_regular, 
                             paste("Regular Expression Distances -", tissue)),
                create_boxplot(df_log2_tissue, type_combinations_log2, "wilcox.test", y_breaks_log2, 
                             paste("Log2 Expression Distances -", tissue))
            )
            
            output_pdf <- file.path(results_dir, paste0(filename_prefix, "_", tissue, "_combined.pdf"))
            ggsave(output_pdf,
                   marrangeGrob(plot_list, nrow=2, ncol=2, top=""),
                   width = 16, height = 12,
                   device = "pdf")
        }
    }
}






process_regular_distances <- function(data_pattern, loaded_objects) {
    data_names <- loaded_objects[grepl(data_pattern, loaded_objects)]
    
    data_names %>%
        map_df(function(data_name) {
            current_data <- get(data_name)
            type_name <- sub(data_pattern, "", data_name)
            
            tibble(
                Type = type_name,
                Distance = unlist(current_data)
            )
        })
}

process_tissue_distances <- function(data_pattern, loaded_objects) {
    data_names <- loaded_objects[grepl(data_pattern, loaded_objects)]
    
    data_names %>%
        map_df(function(data_name) {
            current_data <- get(data_name)
            type_name <- sub(data_pattern, "", data_name)
            
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
}








# Function to process tissue statistics data
process_tissue_statistics <- function(file_path, stats_pattern, type_suffix) {
    load(file_path)
    loaded_objects <- ls()
    data_names_tissue <- loaded_objects[grepl(stats_pattern, loaded_objects)]
    
    df_mean.dists_tissue <- data.frame()
    df_median.dists_tissue <- data.frame()
    
    for (data_name in data_names_tissue) {
        current_data <- get(data_name)
        type_name <- sub(type_suffix, "", data_name)
        
        temp_mean_df <- current_data %>%
            select(Family, Tissue, Mean) %>%
            rename(Distance = Mean) %>%
            mutate(Type = type_name) %>%
            select(Family, Type, Tissue, Distance)
        
        temp_median_df <- current_data %>%
            select(Family, Tissue, Median) %>%
            rename(Distance = Median) %>%
            mutate(Type = type_name) %>%
            select(Family, Type, Tissue, Distance)
        
        df_mean.dists_tissue <- bind_rows(df_mean.dists_tissue, temp_mean_df)
        df_median.dists_tissue <- bind_rows(df_median.dists_tissue, temp_median_df)
    }
    
    return(list(mean = df_mean.dists_tissue, median = df_median.dists_tissue))
}


# Function to create and save tissue boxplots
create_tissue_boxplots <- function(data, metric_type, test_type, is_log2, results_dir) {
    tissue_types <- unique(data$Tissue)
    plot_list <- list()
    
    log2_suffix <- if(is_log2) "log2_" else ""
    
    for (tissue in tissue_types) {
        df_tissue <- subset(data, Tissue == tissue)
        
        boxplot_tissue <- ggplot(df_tissue, aes(x = Type, y = Distance, fill = Type)) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
            labs(
                title = paste(metric_type, "Expression Distances", if(is_log2) "(log2)" else "", 
                            paste0("(", test_type, ") -"), tissue),
                y = "Distance",
                x = "Gene Type"
            ) +
            theme_pubr(border = TRUE) +
            scale_y_continuous(breaks = seq(0, max(df_tissue$Distance), by = 0.1)) +
            theme(
                plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
                axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
                axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
                axis.text.x = element_text(size = 10),
                plot.margin = margin(r = 30)
            )
        
        types <- unique(df_tissue$Type)
        if (length(types) >= 2) {
            type_combinations <- combn(types, 2, simplify = FALSE)
            boxplot_tissue <- boxplot_tissue +
                geom_signif(
                    comparisons = type_combinations,
                    test = test_type,
                    test.args = list(alternative = "two.sided"),
                    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
                    step_increase = 0.05,
                    tip_length = 0.005,
                    color = "black",
                    size = 0.3,
                    textsize = 2.5
                )
        }
        
        plot_list[[tissue]] <- boxplot_tissue
    }
    
    output_pdf <- file.path(results_dir, 
                           paste0("tissues_", tolower(metric_type), "_boxplot_combined_", 
                                 log2_suffix, "(", test_type, ").pdf"))
    ggsave(output_pdf, marrangeGrob(plot_list, nrow=1, ncol=1, top=""), 
           width = 12, height = 8)
}