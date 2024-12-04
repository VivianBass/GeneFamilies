#' Create Statistical Boxplots with Annotations
#'
#' @param data A data frame containing columns 'Type' and 'Distance'
#' @param type_combinations List of vectors, each containing two group names to compare
#' @param test_type Character string specifying the statistical test
#' @param y_breaks Numeric vector specifying the breaks for y-axis
#' @param title_prefix Character string for plot title prefix
#'
#' @return A ggplot object with annotated boxplot
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter stat_summary
#' @importFrom ggpubr theme_pubr
boxplots <- function(data, type_combinations, test_type, y_breaks, title_prefix = "Complete Expression Distances") {
    # Calculate counts per type
    counts <- table(data$Type)
    
    # Calculate the mean of means for each group using the Distance column
    group_means <- tapply(data$Distance, data$Type, mean, na.rm = TRUE)
    hline_value <- mean(group_means, na.rm = TRUE)

    ggplot(data, aes(x = Type, y = Distance, fill = Type)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
        stat_summary(
            fun = mean,
            geom = "text",
            aes(label = sprintf("%.2f", after_stat(y))),
            position = position_dodge(width = 0.75),
            size = 2.5,
            color = "#470050",
            vjust = 0.5,
            hjust = -2.75
        ) +
        stat_summary(
            fun = mean,
            geom = "errorbar",
            aes(ymin = after_stat(y), ymax = after_stat(y)),
            width = 0.5,
            linetype = "dashed",
            linewidth = 1,
            color = "#9b0000"
        ) +
        geom_hline(
            yintercept = hline_value,
            color = "red",
            linetype = "dotted",
            linewidth = 0.5,
            alpha = 0.7
        ) +
        annotate(
            "text",
            x = -Inf,
            y = hline_value,
            label = sprintf("Mean: %.2f", hline_value),
            hjust = -0.05,
            vjust = -0.5,
            color = "#470050",
            size = 2.5
        ) +
        labs(
            title = paste(title_prefix, paste0("(", test_type, ")")), 
            y = "Distance",
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

#' Create Combined Tissue Boxplots for Mean and Median Data
#'
#' @param mean_data Data frame containing mean expression distances
#' @param median_data Data frame containing median expression distances
#' @param log2_mean_data Data frame containing log2 mean distances
#' @param log2_median_data Data frame containing log2 median distances
#' @param results_dir Directory path for saving output files
#'
#' @return None (saves plots to files)
#'
#' @importFrom gridExtra marrangeGrob
#' @importFrom ggplot2 ggsave
create_tissue_boxplots_combined <- function(mean_data, median_data, log2_mean_data, log2_median_data, results_dir) {
    tissue_types <- unique(mean_data$Tissue)
    
    for (tissue in tissue_types) {
        # Filter data for current tissue
        mean_tissue <- subset(mean_data, Tissue == tissue)
        median_tissue <- subset(median_data, Tissue == tissue)
        log2_mean_tissue <- subset(log2_mean_data, Tissue == tissue)
        log2_median_tissue <- subset(log2_median_data, Tissue == tissue)
        
        types <- unique(mean_tissue$Type)
        type_combinations <- combn(types, 2, simplify = FALSE)
        
        # Calculate y-breaks
        y_breaks_mean <- seq(0, max(mean_tissue$Distance, na.rm = TRUE), by = 0.2)
        y_breaks_median <- seq(0, max(median_tissue$Distance, na.rm = TRUE), by = 0.2)
        y_breaks_log2_mean <- seq(0, max(log2_mean_tissue$Distance, na.rm = TRUE), by = 0.2)
        y_breaks_log2_median <- seq(0, max(log2_median_tissue$Distance, na.rm = TRUE), by = 0.2)
        
        # Regular data plots
        regular_plots <- list(
            boxplots(mean_tissue, type_combinations, "t.test", y_breaks_mean,
                    paste("Mean Expression Distances -", tissue)),
            boxplots(mean_tissue, type_combinations, "wilcox.test", y_breaks_mean,
                    paste("Mean Expression Distances -", tissue)),
            boxplots(median_tissue, type_combinations, "t.test", y_breaks_median,
                    paste("Median Expression Distances -", tissue)),
            boxplots(median_tissue, type_combinations, "wilcox.test", y_breaks_median,
                    paste("Median Expression Distances -", tissue))
        )
        
        # Log2 data plots
        log2_plots <- list(
            boxplots(log2_mean_tissue, type_combinations, "t.test", y_breaks_log2_mean,
                    paste("Mean Log2 Expression Distances -", tissue)),
            boxplots(log2_mean_tissue, type_combinations, "wilcox.test", y_breaks_log2_mean,
                    paste("Mean Log2 Expression Distances -", tissue)),
            boxplots(log2_median_tissue, type_combinations, "t.test", y_breaks_log2_median,
                    paste("Median Log2 Expression Distances -", tissue)),
            boxplots(log2_median_tissue, type_combinations, "wilcox.test", y_breaks_log2_median,
                    paste("Median Log2 Expression Distances -", tissue))
        )
        
        # Save regular plots
        ggsave(file.path(results_dir, paste0("tissue_boxplots_regular_", tissue, "_combined.pdf")),
               marrangeGrob(regular_plots, nrow=2, ncol=2, top=""),
               width = 16, height = 12,
               device = "pdf")
        
        # Save log2 plots
        ggsave(file.path(results_dir, paste0("tissue_boxplots_log2_", tissue, "_combined.pdf")),
               marrangeGrob(log2_plots, nrow=2, ncol=2, top=""),
               width = 16, height = 12,
               device = "pdf")
    }
}

#' Create Combined Tissue Boxplots for Regular and Log2 Data
#'
#' @param regular_data Data frame containing regular expression distances
#' @param log2_data Data frame containing log2 expression distances
#' @param results_dir Directory path for saving output files
#' @param filename_prefix Prefix for output filenames
#'
#' @return None (saves plots to files)
#'
#' @importFrom dplyr filter
#' @importFrom gridExtra marrangeGrob
#' @importFrom ggplot2 ggsave
create_tissue_boxplots_all_combined <- function(regular_data, log2_data, results_dir, filename_prefix) {
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
                boxplots(df_regular_tissue, type_combinations_regular, "t.test", y_breaks_regular, 
                             paste("Regular Expression Distances -", tissue)),
                boxplots(df_log2_tissue, type_combinations_log2, "t.test", y_breaks_log2, 
                             paste("Log2 Expression Distances -", tissue)),
                boxplots(df_regular_tissue, type_combinations_regular, "wilcox.test", y_breaks_regular, 
                             paste("Regular Expression Distances -", tissue)),
                boxplots(df_log2_tissue, type_combinations_log2, "wilcox.test", y_breaks_log2, 
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


