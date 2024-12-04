
#' Calculate Angles for a Gene Group
#'
#' This function calculates diagnostic angles for a given set of genes based on RNA-seq expression profiles.
#'
#' @param genes A list of gene groups where each element is a vector of gene identifiers.
#' @param rna.seq.exp.profils A data frame containing RNA-seq expression profiles. 
#'   Must include a column `FBpp_ID` for gene identifiers and additional columns for tissue-specific expression levels.
#' @param tissues A vector of column names representing tissues in `rna.seq.exp.profils` for angle calculation.
#'
#' @return A data frame containing the following columns:
#' \describe{
#'   \item{FBpp_ID}{Gene identifiers matching the input `genes`.}
#'   \item{angle.diag}{Calculated diagnostic angle values for each gene.}
#' }
#' If no matching genes are found, a warning is issued and an empty data frame is returned.
#'
#' @details
#' The function computes angles using the `cosDiag` function applied to tissue expression levels for each gene, normalized by \eqn{\sqrt{2}}.
#' Only genes present in both the `genes` list and the `FBpp_ID` column of `rna.seq.exp.profils` are processed.
#' NA or invalid values are filtered out from the resulting data frame.
#'
#' @note Ensure that the `cosDiag` function is defined in your environment and accepts the appropriate input format.
#'
#' @importFrom dplyr filter
#' @importFrom parallel mclapply
#'
#' @examples
#' \dontrun{
#' # Example data
#' genes <- list(group1 = c("gene1", "gene2"), group2 = c("gene3", "gene4"))
#' rna.seq.exp.profils <- data.frame(
#'   FBpp_ID = c("gene1", "gene2", "gene3"),
#'   tissue1 = c(1.2, 0.8, 1.1),
#'   tissue2 = c(0.9, 1.0, 0.7)
#' )
#' tissues <- c("tissue1", "tissue2")
#'
#' # Compute angles
#' calculate_angles(genes, rna.seq.exp.profils, tissues)
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

#' Create Angle or Versatility Plot
#'
#' Creates a boxplot visualization for either expression angle to diagonal or relative expression versatility,
#' including statistical significance comparisons between gene types.
#'
#' @param data A data frame containing the following columns:
#'   \itemize{
#'     \item gene.type: Factor indicating the type of gene
#'     \item angle.diag or rel.vers: Numeric values for plotting (depending on plot_type)
#'   }
#' @param plot_type Character string specifying the type of plot to create:
#'   \itemize{
#'     \item "angle": Creates an Expression Angle to Diagonal plot
#'     \item "versatility": Creates a Relative Expression Versatility plot
#'   }
#' @param test_type Character string specifying the statistical test to use (e.g., "wilcox.test")
#'
#' @return A ggplot object containing the generated plot
#'
#' @details
#' The function creates a boxplot with jittered points and statistical significance indicators.
#' For angle plots, the statistical test uses a "greater" alternative hypothesis,
#' while versatility plots use a "two.sided" alternative.
#'
#' @note
#' Requires the following global variables to be defined:
#' \itemize{
#'   \item type_combinations: List of gene type pairs for statistical comparison
#'   \item results_dir: Directory path for saving the plot
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter labs scale_y_continuous theme ggsave
#' @importFrom ggpubr theme_pubr
#' @importFrom ggsignif geom_signif
#'
#' @examples
#' \dontrun{    
#' data <- data.frame(
#'   gene.type = factor(rep(c("TypeA", "TypeB"), each = 100)),
#'   angle.diag = rnorm(200),
#'   rel.vers = rnorm(200)
#' )
#' type_combinations <- list(c("TypeA", "TypeB"))
#' results_dir <- "path/to/results"
#' 
#' # Create angle plot
#' plot1 <- create_angle_versatility_plot(data, "angle", "wilcox.test")
#' 
#' # Create versatility plot
#' plot2 <- create_angle_versatility_plot(data, "versatility", "wilcox.test")
#' }
#'
#' @export
create_angle_versatility_plot <- function(data, plot_type, test_type) {
    # Set up variables based on plot type
    if (plot_type == "angle") {
        y_var <- "angle.diag"
        title <- "Expression Angle To Diagonal"
        y_label <- "Relative Tissue Specificity"
        output_name <- "expressionAngleToDiagonalBoxplot"
    } else {
        y_var <- "rel.vers"
        title <- "Relative Expression Versatility"
        y_label <- "Relative Tissue Versatility"
        output_name <- "relativeExpressionVersatilityBoxplot"
    }
    
    # Create plot
    plot <- ggplot(data, aes(x = gene.type, y = .data[[y_var]], fill = gene.type)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
        labs(title = title, y = y_label, x = "Gene Type") +
        theme_pubr(border = TRUE) +
        scale_y_continuous(breaks = seq(0, max(data[[y_var]], na.rm = TRUE), by = 0.1)) +
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
    
    # Save plot
    ggsave(file.path(results_dir, paste0(output_name, "_", test_type, ".pdf")),
           plot, width = 10, height = 8)
    
    return(plot)
}



#' Create Angle or Versatility Plot
#'
#' Creates a boxplot visualization for either expression angle to diagonal or relative expression versatility,
#' including statistical significance comparisons between gene types.
#'
#' @param data A data frame containing the following columns:
#'   \itemize{
#'     \item gene.type: Factor indicating the type of gene
#'     \item angle.diag or rel.vers: Numeric values for plotting (depending on plot_type)
#'   }
#' @param plot_type Character string specifying the type of plot to create:
#'   \itemize{
#'     \item "angle": Creates an Expression Angle to Diagonal plot
#'     \item "versatility": Creates a Relative Expression Versatility plot
#'   }
#' @param test_type Character string specifying the statistical test to use (e.g., "wilcox.test")
#'
#' @return A ggplot object containing the generated plot
#'
#' @details
#' The function creates a boxplot with jittered points and statistical significance indicators.
#' For angle plots, the statistical test uses a "greater" alternative hypothesis,
#' while versatility plots use a "two.sided" alternative.
#'
#' @note
#' Requires the following global variables to be defined:
#' \itemize{
#'   \item type_combinations: List of gene type pairs for statistical comparison
#'   \item results_dir: Directory path for saving the plot
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter labs scale_y_continuous theme ggsave
#' @importFrom ggpubr theme_pubr
#' @importFrom ggsignif geom_signif
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   gene.type = factor(rep(c("TypeA", "TypeB"), each = 100)),
#'   angle.diag = rnorm(200),
#'   rel.vers = rnorm(200)
#' )
#' type_combinations <- list(c("TypeA", "TypeB"))
#' results_dir <- "path/to/results"
#' 
#' # Create angle plot
#' plot1 <- create_angle_versatility_plot(data, "angle", "wilcox.test")
#' 
#' # Create versatility plot
#' plot2 <- create_angle_versatility_plot(data, "versatility", "wilcox.test")
#' }
#'
#' @export
