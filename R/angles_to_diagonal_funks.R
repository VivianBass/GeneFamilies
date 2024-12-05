
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
#' @param suffix Suffix to add to the output file names
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
create_angle_versatility_plots <- function(data, test_types, results_dir, suffix = "") {
    angle_plots <- list()
    versatility_plots <- list()
    
    # Setup type combinations for significance testing
    gene_types <- levels(data$gene.type)
    type_combinations <- combn(gene_types, 2, simplify = FALSE)
    
    for (test in test_types) {
        # Assign type_combinations to global environment for create_single_plot
        assign("type_combinations", type_combinations, envir = .GlobalEnv)
        
        angle_plots[[test]] <- create_single_plot(data, "angle", test)
        versatility_plots[[test]] <- create_single_plot(data, "versatility", test)
    }
    
    # Combine and save angle plots
    combined_angle_plots <- grid.arrange(
        grobs = list(angle_plots[["t.test"]], angle_plots[["wilcox.test"]]),
        ncol = 2,
        top = paste0("Expression Angle Plots", suffix)
    )
    ggsave(
        file.path(results_dir, paste0("combined_angle_plots", suffix, ".pdf")),
        combined_angle_plots,
        width = 20, height = 8
    )
    
    # Combine and save versatility plots
    combined_versatility_plots <- grid.arrange(
        grobs = list(versatility_plots[["t.test"]], versatility_plots[["wilcox.test"]]),
        ncol = 2,
        top = paste0("Expression Versatility Plots", suffix)
    )
    ggsave(
        file.path(results_dir, paste0("combined_versatility_plots", suffix, ".pdf")),
        combined_versatility_plots,
        width = 20, height = 8
    )
}


# --------------------------------------------------------------------------------

#' Perform Statistical Tests on Gene Expression Data
#' 
#' FUNCTION FOR ONE DIRECTIONAL TESTS
#'
#' @param data Data frame containing gene expression data with columns:
#'   - gene.type: Factor indicating gene groups for comparison
#'   - Additional columns specified by column_name parameter
#' @param valid_groups List containing valid gene groups for comparison
#' @param column_name String specifying which column to analyze (e.g., "angle.diag" or "rel.vers")
#'
#' @return List containing two elements:
#'   - t_test: Results of t-test with adjusted p-values and significance levels
#'   - wilcox: Results of Wilcoxon test with adjusted p-values and significance levels
#'   Returns NULL if fewer than 2 valid groups are available
#'
#' @details
#' Performs both t-test and Wilcoxon test with "greater" alternative hypothesis.
#' P-values are adjusted using Benjamini-Hochberg method.
#' Significance levels are coded as:
#'   - "ns": p >= 0.05
#'   - "*": p < 0.05
#'   - "**": p < 0.01
#'   - "***": p < 0.001
#'
#' @importFrom dplyr filter mutate case_when
#' @importFrom rstatix t_test wilcox_test adjust_pvalue
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   gene.type = factor(c("type1", "type2", "type1", "type2")),
#'   angle.diag = c(0.5, 0.7, 0.6, 0.8)
#' )
#' valid_groups <- list(groups = c("type1", "type2"))
#' results <- perform_tests(data, valid_groups, "angle.diag")
#' }
perform_tests_X <- function(data, valid_groups, column_name) {
    if (length(valid_groups[[1]]) >= 2) {
        formula <- as.formula(paste(column_name, "~ gene.type"))
        
        t_test_result <- data %>%
            filter(gene.type %in% valid_groups[[1]]) %>%
            t_test(formula, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(analysis = column_name, 
                   test_type = "t-test",
                   p.adj.signif = case_when(
                       p.adj >= 0.05 ~ "ns",
                       p.adj < 0.001 ~ "***",
                       p.adj < 0.01 ~ "**",
                       p.adj < 0.05 ~ "*"
                   ))
        
        wilcox_result <- data %>%
            filter(gene.type %in% valid_groups[[1]]) %>%
            wilcox_test(formula, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(analysis = column_name, 
                   test_type = "wilcox",
                   p.adj.signif = case_when(
                       p.adj >= 0.05 ~ "ns",
                       p.adj < 0.001 ~ "***",
                       p.adj < 0.01 ~ "**",
                       p.adj < 0.05 ~ "*"
                   ))
        
        return(list(t_test = t_test_result, wilcox = wilcox_result))
    }
    return(NULL)
}

#' Perform Statistical Tests on Gene Expression Data
#'
#' Conducts t-tests and Wilcoxon tests on all possible pairs of groups in both directions
#' (e.g., group1 vs group2 and group2 vs group1) for a specified metric column.
#'
#' @param data Data frame containing 'gene.type' and the metric column to be tested.
#' @param valid_groups List of valid groups for analysis, each as a character vector.
#' @param column_name Name of the metric column to perform tests on.
#'
#' @return A list with t-test and Wilcoxon test results for each group pair, or NULL if no valid results.
#'
#' @details
#' - Requires at least 2 observations per group.
#' - Adjusts p-values using the Benjamini-Hochberg method.
#' - Tests each group pair in both directions, ensuring the order of comparison is considered.
#'
#' @importFrom dplyr filter mutate arrange bind_rows
#' @importFrom rstatix t_test wilcox_test adjust_pvalue
#'
#' @examples
#' data <- data.frame(gene.type = c("A", "A", "B", "B"), Distance = c(1.2, 1.5, 2.3, 2.1))
#' results <- perform_tests(data, list(c("A", "B")), "Distance")
#'
#' @export
perform_tests <- function(data, valid_groups, column_name) {
    if (length(valid_groups[[1]]) >= 2) {
        # Create all possible combinations of groups in both directions
        groups <- valid_groups[[1]]
        group_pairs <- expand.grid(group1 = groups, group2 = groups)
        # Remove self-comparisons
        group_pairs <- group_pairs[group_pairs$group1 != group_pairs$group2,]
        
        # Initialize empty results
        t_test_results <- list()
        wilcox_results <- list()
        
        formula <- as.formula(paste(column_name, "~ gene.type"))
        
        # Loop through each pair
        for(i in 1:nrow(group_pairs)) {
            g1 <- group_pairs$group1[i]
            g2 <- group_pairs$group2[i]
            
            # Get data for current pair and ensure correct order for comparison
            test_data <- data %>% 
                filter(gene.type %in% c(g1, g2)) %>%
                mutate(gene.type = factor(gene.type, levels = c(g2, g1)))  # Order matters for alternative="greater"
            
            # T-test
            t_result <- test_data %>%
                t_test(formula, alternative = "greater") %>%
                mutate(
                    analysis = column_name,
                    test_type = "t-test",
                    group1 = g1,
                    group2 = g2
                )
            t_test_results[[i]] <- t_result
            
            # Wilcoxon test
            w_result <- test_data %>%
                wilcox_test(formula, alternative = "greater") %>%
                mutate(
                    analysis = column_name,
                    test_type = "wilcox",
                    group1 = g1,
                    group2 = g2
                )
            wilcox_results[[i]] <- w_result
        }
        
        # Combine results
        t_test_result <- bind_rows(t_test_results) %>%
            adjust_pvalue(method = "BH") %>%
            mutate(
                p.adj.signif = case_when(
                    p.adj >= 0.05 ~ "ns",
                    p.adj < 0.001 ~ "***",
                    p.adj < 0.01 ~ "**",
                    p.adj < 0.05 ~ "*"
                )
            ) %>%
            arrange(test_type)
        
        wilcox_result <- bind_rows(wilcox_results) %>%
            adjust_pvalue(method = "BH") %>%
            mutate(
                p.adj.signif = case_when(
                    p.adj >= 0.05 ~ "ns",
                    p.adj < 0.001 ~ "***",
                    p.adj < 0.01 ~ "**",
                    p.adj < 0.05 ~ "*"
                )
            ) %>%
            arrange(test_type)
        
        return(list(t_test = t_test_result, wilcox = wilcox_result))
    }
    return(NULL)
}


#' Perform Statistical Tests for Two Data Frames
#'
#' @param df1 First data frame with columns 'gene.type', 'angle.diag', and 'rel.vers'
#' @param df2 Second data frame with same structure as df1
#' @param output_file Name of CSV file to save test results
#'
#' @return List containing test results:
#'   - angle_t_test: T-test results for angle.diag
#'   - angle_wilcox: Wilcoxon test results for angle.diag
#'   - relvers_t_test: T-test results for rel.vers
#'   - relvers_wilcox: Wilcoxon test results for rel.vers
#'   - summary: Combined test results
#'
#' @importFrom dplyr bind_rows group_by summarise filter pull
#' @importFrom utils write.csv
#'
#' @examples
#' \dontrun{
#' results <- perform_statistical_tests_for_columns(plot.df, plot.df_log2, "test_results.csv")
#' }
perform_statistical_tests_for_columns <- function(df1, df2, output_file) {
    # Combine data frames and add a source column
    combined_df <- bind_rows(
        df1 %>% mutate(source = "df1"),
        df2 %>% mutate(source = "df2")
    )
    
    # Check for valid groups
    valid_groups <- combined_df %>%
        group_by(gene.type) %>%
        summarise(n = n(), .groups = 'drop') %>%
        filter(n > 1) %>%
        pull(gene.type)
    
    # Perform tests
    test_results <- tryCatch({
        # Perform tests for angle.diag
        angle_tests <- perform_tests(combined_df, list(angle = valid_groups), "angle.diag")
        
        # Perform tests for rel.vers
        relvers_tests <- perform_tests(combined_df, list(relvers = valid_groups), "rel.vers")
        
        if (!is.null(angle_tests) && !is.null(relvers_tests)) {
            test_summary <- bind_rows(
                angle_tests$t_test,
                angle_tests$wilcox,
                relvers_tests$t_test,
                relvers_tests$wilcox
            )
            
            write.csv(test_summary,
                     file.path(results_dir, output_file),
                     row.names = FALSE)
            message("Statistical tests summary exported to CSV")
        }
        
        list(
            angle_t_test = if(!is.null(angle_tests)) angle_tests$t_test else NULL,
            angle_wilcox = if(!is.null(angle_tests)) angle_tests$wilcox else NULL,
            relvers_t_test = if(!is.null(relvers_tests)) relvers_tests$t_test else NULL,
            relvers_wilcox = if(!is.null(relvers_tests)) relvers_tests$wilcox else NULL,
            summary = if(exists("test_summary")) test_summary else NULL
        )
    }, error = function(e) {
        message("Error in statistical tests: ", e$message)
        return(NULL)
    })
    
    return(test_results)
}

# --------------------------------------------------------------------------------

#' Convert Radians to Degrees
#'
#' This function converts angle measurements from radians to degrees.
#'
#' @param radians A numeric value or vector of values in radians
#' @return A numeric value or vector of values in degrees
#' @examples
#' rad2deg(pi)  # returns 180
#' rad2deg(c(pi/2, pi/4))  # returns c(90, 45)
#' @export
rad2deg <- function(radians) {
  return(radians * (180 / pi))
}

#' Convert Radians to Degrees in Nested Lists
#'
#' This function recursively converts all numeric values in a nested list structure
#' from radians to degrees. Non-numeric values are left unchanged.
#'
#' @param data A nested list containing numeric values in radians
#' @return A nested list with the same structure where all numeric values
#'         have been converted from radians to degrees
#' @examples
#' data <- list(
#'   a = list(x = pi, y = "NA"),
#'   b = list(z = pi/2)
#' )
#' convert_radians_to_degrees(data)
#' @export
convert_radians_to_degrees <- function(data) {
  lapply(data, function(sublist) {
    lapply(sublist, function(value) {
      if (is.numeric(value)) {
        return(rad2deg(value))
      } else {
        return(value)
      }
    })
  })
}