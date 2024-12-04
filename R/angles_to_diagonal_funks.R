
#' Calculate Diagnostic Angles for Gene Groups
#'
#' @param genes List of gene groups with gene identifiers
#' @param rna.seq.exp.profils Data frame with RNA-seq profiles (requires FBpp_ID column)
#' @param tissues Vector of tissue column names for angle calculation
#'
#' @return Data frame with columns:
#'   - FBpp_ID: Gene identifiers
#'   - angle.diag: Calculated diagnostic angles
#'
#' @details Calculates angles using cosDiag function normalized by sqrt(2)
#'
#' @importFrom dplyr filter
#' @importFrom parallel mclapply
#'
#' @examples
#' \dontrun{
#' genes <- list(g1 = c("gene1", "gene2"))
#' profiles <- data.frame(
#'   FBpp_ID = c("gene1", "gene2"),
#'   t1 = c(1.2, 0.8),
#'   t2 = c(0.9, 1.0)
#' )
#' angles <- calculate_angles(genes, profiles, c("t1", "t2"))
#' }
#'
#' @export
calculate_angles <- function(genes, rna.seq.exp.profils, tissues) {
  genes.expr <- intersect(unlist(genes), rna.seq.exp.profils$FBpp_ID)
  
  if(length(genes.expr) == 0) {
    warning("No matching genes found")
    return(data.frame())
  }

  expr.angle.diag.df <- data.frame(
    FBpp_ID = genes.expr,
    angle.diag = as.numeric(mclapply(genes.expr, function(x) {
      cosDiag(rna.seq.exp.profils[which(rna.seq.exp.profils$FBpp_ID == x), tissues])/sqrt(2)
    })),
    stringsAsFactors = FALSE
  )

  expr.angle.diag.df %>%
    filter(!is.na(angle.diag) & angle.diag != "NA" & angle.diag != "")
}

#' Validate Data Frames and Retain Original Names
#'
#' This function validates a list of data frames, keeping only those that are non-null and have at least one row. 
#' The names of the valid data frames in the input list are preserved in the returned list.
#'
#' @param df_list A named list of data frames to be validated.
#'
#' @return A named list of validated data frames. 
#'   Only data frames that are non-null and contain at least one row are included in the output list.
#'   If a data frame is null or empty, a warning message is displayed.
#'
#' @details
#' This function iterates over the input list of data frames, checks each data frame for validity, 
#' and retains the original name of each valid data frame in the output list.
#'
#' @examples
#' # Example usage
#' df_list <- list(
#'   df1 = data.frame(a = 1:3, b = 4:6),
#'   df2 = NULL,
#'   df3 = data.frame()
#' )
#'
#' validate_angle_dataframes(df_list)
#' # Output: List containing only 'df1', with a warning for 'df2' and 'df3'.
#'
#' @export
validate_angle_dataframes <- function(df_list) {
    valid_dfs <- list()
    
    for (name in names(df_list)) {
        df <- df_list[[name]]
        if (!is.null(df) && nrow(df) > 0) {
            valid_dfs[[name]] <- df
        } else {
            message(sprintf("Warning: %s dataframe is empty or null", name))
        }
    }
    
    if (length(valid_dfs) == 0) {
        return(structure(list(), class = "list"))
    }
    
    return(valid_dfs)
}

#' Create Boxplot for Expression Angle or Versatility
#'
#' Generates a boxplot for expression angle or versatility with significance annotations.
#'
#' @param data Data frame with 'gene.type', 'angle.diag', and 'rel.vers' columns.
#' @param plot_type "angle" for expression angle or "versatility" for relative versatility.
#' @param test_type Statistical test: "t.test" or "wilcox.test".
#'
#' @return A ggplot object with the boxplot and significance annotations.
#'
#' @note Requires global variable 'type_combinations' for comparisons.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter labs theme
#' @importFrom ggpubr theme_pubr
#' @importFrom ggsignif geom_signif
#'
#' @examples
#' \dontrun{
#' data <- data.frame(gene.type = factor(rep(c("A", "B"), each = 100)), angle.diag = rnorm(200))
#' plot <- create_single_plot(data, "angle", "wilcox.test")
#' }
#'
#' @export
create_single_plot <- function(data, plot_type, test_type) {
    # Debug: Print column names
    print("Column names in data:")
    print(names(data))
    
    # Calculate counts per type
    counts <- table(data$gene.type)
    
    # Set up variables based on plot type and test type
    if (plot_type == "angle") {
        y_var <- "angle.diag"
        title <- sprintf("Expression Angle To Diagonal (%s)", 
                        ifelse(test_type == "t.test", "T-test", "Wilcoxon"))
        y_label <- "Relative Tissue Specificity"
    } else {
        y_var <- "rel.vers"
        title <- sprintf("Relative Expression Versatility (%s)", 
                        ifelse(test_type == "t.test", "T-test", "Wilcoxon"))
        y_label <- "Relative Tissue Versatility"
    }

    # Calculate the mean of means for each group
    group_means <- tapply(data[[y_var]], data$gene.type, mean, na.rm = TRUE)
    hline_value <- mean(group_means, na.rm = TRUE)
    
    # Create plot
    plot <- ggplot(data, aes(x = gene.type, y = .data[[y_var]], fill = gene.type)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
        stat_summary(
            fun = mean,
            geom = "text",
            aes(label = sprintf("%.2f", after_stat(y))),
            position = position_dodge(width = 0.75),
            size = 3,
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
            color = "#470050",
            size = 3
        ) +
        labs(title = title, y = y_label, x = "Gene Type") +
        theme_pubr(border = TRUE) +
        scale_y_continuous(breaks = seq(0, max(data[[y_var]], na.rm = TRUE), by = 0.1)) +
        scale_x_discrete(labels = function(x) paste0(x, "\n(n=", counts[x], ")")) +
        theme(
            plot.title = element_text(size = 12, face = "bold", margin = margin(t = 20, b = 20), hjust = 0.5),
            axis.title.x = element_text(size = 10, margin = margin(t = 20, b = 20), hjust = 0.5),
            axis.title.y = element_text(size = 10, margin = margin(t = 20, r = 20, b = 20, l = 20)),
            axis.text.x = element_text(size = 10),
            plot.margin = margin(r = 30)
        ) +
        geom_signif(
            comparisons = type_combinations,
            test = test_type,
            test.args = list(alternative = if(plot_type == "angle" || test_type == "wilcox.test") "greater" else "two.sided"),
            map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
            step_increase = 0.05,
            tip_length = 0.005,
            vjust = 0.5,
            color = "black",
            size = 0.3,
            textsize = 2.5
        )
    
    return(plot)
}

#' Create Combined Angle and Versatility Plot Sets
#'
#' Creates and saves two combined plot sets: one for expression angles and one for versatility,
#' each combining t-test and Wilcoxon test results side by side.
#'
#' @param data Data frame with gene expression data
#' @param test_types Vector of statistical tests to perform
#' @param results_dir Directory path for saving output PDFs
#'
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggsave
#'
#' @examples
#' \dontrun{
#' create_angle_versatility_plots(data, c("t.test", "wilcox.test"), "path/to/results")
#' }
#'
#' @export
create_angle_versatility_plots <- function(data, test_types, results_dir) {
    angle_plots <- list()
    versatility_plots <- list()
    
    for (test in test_types) {
        angle_plots[[test]] <- create_single_plot(data, "angle", test)
        versatility_plots[[test]] <- create_single_plot(data, "versatility", test)
    }
    
    # Combine and save angle plots
    combined_angle_plots <- grid.arrange(
        grobs = list(angle_plots[["t.test"]], angle_plots[["wilcox.test"]]),
        ncol = 2,
        top = "Expression Angle Plots"
    )
    ggsave(
        file.path(results_dir, "combined_angle_plots.pdf"),
        combined_angle_plots,
        width = 20, height = 8
    )
    
    # Combine and save versatility plots
    combined_versatility_plots <- grid.arrange(
        grobs = list(versatility_plots[["t.test"]], versatility_plots[["wilcox.test"]]),
        ncol = 2,
        top = "Expression Versatility Plots"
    )
    ggsave(
        file.path(results_dir, "combined_versatility_plots.pdf"),
        combined_versatility_plots,
        width = 20, height = 8
    )
}

