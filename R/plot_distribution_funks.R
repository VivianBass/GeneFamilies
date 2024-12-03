
#' Create Statistical Boxplot for Expression Data
#'
#' @param data Data frame with Type column and metric values
#' @param type_combinations List of type pairs for comparison
#' @param metric Character: "Distance" or "Angle"
#' @param title_prefix Character: Plot title prefix
#' @param test_type Character: Statistical test ("t.test" or "wilcox.test")
#' @param alternative Character: Test direction ("greater", "less", "two.sided")
#' @param y_breaks Numeric vector: Y-axis breaks (optional)
#'
#' @return ggplot object with boxplot and statistical annotations
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter stat_summary
#' @importFrom ggpubr theme_pubr
#' @importFrom ggsignif geom_signif
#'
#' @examples
#' \dontrun{
#' plot <- distance_boxplot(
#'   data = data.frame(Type = c("A", "B"), Distance = c(1, 2)),
#'   type_combinations = list(c("A", "B"))
#' )
#' }
#'
#' @export
distance_boxplot <- function(data,
                             type_combinations,
                             metric = "Distance",
                             title_prefix = "Expression Distances",
                             test_type = "t.test",
                             alternative = "greater",
                             y_breaks = NULL) {
    
    # Rename Distance column to Angle if needed
    if (metric == "Angle" && "Distance" %in% names(data)) {
        data$Angle <- data$Distance
    }
    
    # Validate data has valid values for the metric
    if (all(is.na(data[[metric]]))) {
        stop(paste("No valid values found in", metric, "column"))
    }
    
    # Get valid range for y-axis
    valid_values <- data[[metric]][!is.na(data[[metric]])]
    if (length(valid_values) == 0) {
        y_breaks <- seq(0, 1, 0.2)  # Default range if no valid values
    } else {
        y_breaks <- seq(0, max(valid_values), by = 0.2)
    }
    
    # Calculate counts per type
    counts <- table(data$Type)
    
    # Calculate y_breaks if not provided
    if (is.null(y_breaks)) {
        y_breaks <- seq(0, max(data[[metric]], na.rm = TRUE), by = 0.2)
    }

    hline_value <- mean(aggregate(as.formula(paste(metric, "~ Type")), data = data, FUN = mean)[[metric]])
    
    ggplot(data, aes(x = Type, y = .data[[metric]], fill = Type)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
        stat_summary(
            fun = mean,
            geom = "text",
            aes(label = sprintf("%.2f", after_stat(y))),
            position = position_nudge(x = 0.5, y = 0),
            size = 3,
            color = "red"
        ) +
        geom_hline(
            yintercept = hline_value,
            color = "red",
            linetype = "solid",
            linewidth = 0.5
        ) +
        annotate(
            "text",
            x = -Inf,
            y = hline_value,
            label = sprintf("Mean: %.2f", hline_value),
            hjust = -0.05,
            vjust = -0.5,
            color = "#009c22",
            size = 3
        ) +
        labs(
            title = paste(title_prefix, paste0("(", test_type, ")")),
            y = metric,
            fill = "Type: "
        ) +
        theme_pubr(border = TRUE) +
        scale_y_continuous(breaks = y_breaks) +
        scale_x_discrete(labels = function(x) paste0(x, "\n(n=", counts[x], ")")) +
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
            legend.text = element_text(size = 8),
            plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
            axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
            axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
            axis.text.x = element_text(size = 8),
            plot.margin = margin(r = 30)
        )
}




#' Create Multiple Distance/Angle Boxplots with Statistical Tests
#'
#' @param mean_data Data frame containing mean expression data
#' @param median_data Data frame containing median expression data
#' @param type_combinations List of type pairs for statistical comparison
#' @param metric Character string: "Distance" or "Angle"
#' @param alternative Character string for statistical test direction (default: "greater")
#' @param results_dir Directory path for saving output
#' @param filename Output filename (default: "boxplots_expression_distances_all_tests.pdf")
#' @param plot_width Numeric width of output plot (default: 16)
#' @param plot_height Numeric height of output plot (default: 12)
#' @param nrow Number of rows in combined plot (default: 2)
#' @param ncol Number of columns in combined plot (default: 2)
#'
#' @return List of four ggplot objects: mean/median t-tests and Wilcoxon tests
#'
#' @importFrom ggplot2 ggsave
#' @importFrom gridExtra marrangeGrob
#'
#' @examples
#' \dontrun{
#' plots <- create_distance_boxplots(
#'   mean_data = mean_df,
#'   median_data = median_df,
#'   type_combinations = list(c("TypeA", "TypeB")),
#'   results_dir = "output"
#' )
#' }
#'
#' @export
create_distance_boxplots <- function(mean_data,
                                   median_data,
                                   type_combinations,
                                   metric = "Distance",
                                   alternative = "greater",
                                   results_dir,
                                   filename = "boxplots_expression_distances_all_tests.pdf",
                                   plot_width = 16,
                                   plot_height = 12,
                                   nrow = 2,
                                   ncol = 2) {
    
    # Add Angle column if needed
    if (metric == "Angle") {
        mean_data$Angle <- mean_data$Distance
        median_data$Angle <- median_data$Distance
    }
    
    # T-test plots
    boxplot_mean_ttest <- distance_boxplot(
        data = mean_data,
        type_combinations = type_combinations,
        metric = metric,  # Pass metric
        title_prefix = paste("Mean Expression", metric),
        test_type = "t.test",
        alternative = alternative
    )
    
    boxplot_median_ttest <- distance_boxplot(
        data = median_data,
        type_combinations = type_combinations,
        metric = metric,
        title_prefix = paste("Median Expression", metric),
        test_type = "t.test",
        alternative = alternative
    )
    
    boxplot_mean_wilcox <- distance_boxplot(
        data = mean_data,
        type_combinations = type_combinations,
        metric = metric,
        title_prefix = paste("Mean Expression", metric),
        test_type = "wilcox.test",
        alternative = alternative
    )
    
    boxplot_median_wilcox <- distance_boxplot(
        data = median_data,
        type_combinations = type_combinations,
        metric = metric,
        title_prefix = paste("Median Expression", metric),
        test_type = "wilcox.test",
        alternative = alternative
    )
    
    output_pdf <- file.path(results_dir, filename)
    plot_list <- list(boxplot_median_ttest, boxplot_mean_ttest,
                     boxplot_median_wilcox, boxplot_mean_wilcox)
    ggsave(output_pdf,
           marrangeGrob(plot_list, nrow=nrow, ncol=ncol, top=""),
           width = plot_width, height = plot_height,
           device = "pdf")
    
    return(plot_list)
}


# -------------------------------------------------------------------


#' Create Comprehensive Boxplot with Statistical Analysis
#'
#' @param data Data frame containing Type column and metric values
#' @param type_combinations List of type pairs for statistical comparison
#' @param test_type Character string specifying statistical test type
#' @param y_breaks Numeric vector for y-axis breaks
#' @param metric Character string specifying metric column (default: "Angle")
#' @param title_prefix Character string for plot title prefix
#'
#' @return A ggplot object with boxplot, jittered points, mean values, and statistical comparisons
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter stat_summary geom_hline annotate
#' @importFrom ggpubr theme_pubr
#' @importFrom ggsignif geom_signif
#'
#' @examples
#' \dontrun{
#' plot <- distance_boxplots_all(
#'   data = data.frame(Type = c("A", "B"), Angle = c(1, 2)),
#'   type_combinations = list(c("A", "B")),
#'   test_type = "t.test",
#'   y_breaks = seq(0, 2, 0.5)
#' )
#' }
#'
#' @export
distance_boxplots_all <- function(data, 
                                type_combinations, 
                                test_type, 
                                y_breaks, 
                                metric = "Angle",
                                title_prefix = paste("Complete Expression", metric)) {
    
    # Calculate counts per type
    counts <- table(data$Type)
    
    # Use metric parameter for y-axis mapping
    y_col <- sym(metric)

    hline_value <- mean(aggregate(as.formula(paste(metric, "~ Type")), data = data, FUN = mean)[[metric]])
    
    ggplot(data, aes(x = Type, y = !!y_col, fill = Type)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
        stat_summary(
            fun = mean,
            geom = "text",
            aes(label = sprintf("%.2f", after_stat(y))),
            position = position_nudge(x = 0.5, y = 0),
            size = 3,
            color = "red"
        ) +
        geom_hline(
            yintercept = hline_value,
            color = "red",
            linetype = "solid",
            linewidth = 0.5
        ) +
        annotate(
            "text",
            x = -Inf,
            y = hline_value,
            label = sprintf("Mean: %.2f", hline_value),
            hjust = -0.05,
            vjust = -0.5,
            color = "#009c22",
            size = 3
        ) +
        labs(
            title = paste(title_prefix, paste0("(", test_type, ")")),
            y = metric,
            fill = "Type: "
        ) +
        theme_pubr(border = TRUE) +
        scale_y_continuous(breaks = y_breaks) +
        scale_x_discrete(labels = function(x) paste0(x, "\n(n=", counts[x], ")")) +
        geom_signif(
            comparisons = type_combinations,
            test = test_type,
            test.args = list(alternative = "greater"),
            map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
            step_increase = 0.05,
            tip_length = 0.005,
            color = "black",
            size = 0.3,
            textsize = 2.5
        ) +
        theme(
            legend.text = element_text(size = 8),
            plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
            axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
            axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
            axis.text.x = element_text(size = 8),
            plot.margin = margin(r = 30)
        )
}



#' Create Comprehensive Boxplot with Statistical Analysis
#'
#' @param data Data frame containing Type column and metric values
#' @param type_combinations List of type pairs for statistical comparison
#' @param test_type Character string specifying statistical test type
#' @param y_breaks Numeric vector for y-axis breaks
#' @param metric Character string specifying metric column (default: "Angle")
#' @param title_prefix Character string for plot title prefix
#'
#' @return A ggplot object with boxplot, jittered points, mean values, and statistical comparisons
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter stat_summary geom_hline annotate
#' @importFrom ggpubr theme_pubr
#' @importFrom ggsignif geom_signif
#'
#' @examples
#' \dontrun{
#' plot <- distance_boxplots_all(
#'   data = data.frame(Type = c("A", "B"), Angle = c(1, 2)),
#'   type_combinations = list(c("A", "B")),
#'   test_type = "t.test",
#'   y_breaks = seq(0, 2, 0.5)
#' )
#' }
#'
#' @export
create_distance_boxplots_all <- function(regular_data, 
                                         log2_data, 
                                         results_dir, 
                                         filename_prefix,
                                         metric = "Angle") {
    
    # Process regular data
    df_regular_filtered <- regular_data %>% filter(is.finite(!!sym(metric)))
    types_regular <- unique(df_regular_filtered$Type)
    type_combinations_regular <- combn(types_regular, 2, simplify = FALSE)
    y_breaks_regular <- seq(0, max(df_regular_filtered[[metric]], na.rm = TRUE), by = 0.2)
    
    # Process log2 data
    df_log2_filtered <- log2_data %>% filter(is.finite(!!sym(metric)))
    types_log2 <- unique(df_log2_filtered$Type)
    type_combinations_log2 <- combn(types_log2, 2, simplify = FALSE)
    y_breaks_log2 <- seq(0, max(df_log2_filtered[[metric]], na.rm = TRUE), by = 0.2)
    
    # Create all plots
    plot_list <- list(
        distance_boxplots_all(df_regular_filtered, type_combinations_regular, "t.test", y_breaks_regular, metric, paste("Regular Expression", metric)),
        distance_boxplots_all(df_log2_filtered, type_combinations_log2, "t.test", y_breaks_log2, metric, paste("Log2 Expression", metric)),
        distance_boxplots_all(df_regular_filtered, type_combinations_regular, "wilcox.test", y_breaks_regular, metric, paste("Regular Expression", metric)),
        distance_boxplots_all(df_log2_filtered, type_combinations_log2, "wilcox.test", y_breaks_log2, metric, paste("Log2 Expression", metric))
    )
    
    # Save combined plots
    output_pdf <- file.path(results_dir, paste0(filename_prefix, "_combined.pdf"))
    ggsave(output_pdf,
           marrangeGrob(plot_list, nrow=2, ncol=2, top=""),
           width = 16, height = 12,
           device = "pdf")
    
    return(plot_list)
}
