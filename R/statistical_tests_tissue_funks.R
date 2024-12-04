
#' Perform Tissue-Specific Statistical Tests
#'
#' Conducts t-tests and Wilcoxon tests for comparing gene distances between types within tissues. 
#' Results are adjusted for multiple comparisons using the Benjamini-Hochberg method.
#'
#' @param data A dataframe containing distance data with columns `Tissue`, `Type`, and `Distance`.
#' @param valid_groups A named list of dataframes with valid `Tissue` and `Type` combinations for analysis.
#' @param analysis_type A string indicating the analysis type ("mean" or "median").
#' @return A list containing t-test and Wilcoxon test results as dataframes, or `NULL` if no valid groups are available.
#' @examples
#' perform_tissue_tests(data, valid_groups, "mean")
perform_tissue_tests <- function(data, valid_groups, analysis_type) {
    if (nrow(valid_groups[[analysis_type]]) >= 2) {

        t_test_result <- data %>%
            semi_join(valid_groups[[analysis_type]], by = c("Tissue", "Type")) %>%
            group_by(Tissue) %>%
            filter(!is.na(Distance)) %>%
            filter(n_distinct(Type) >= 2) %>%
            t_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(analysis = analysis_type,
                   test_type = "t-test",
                   p.adj.signif = case_when(
                       p >= 0.05 ~ "ns",
                       p < 0.001 ~ "***",
                       p < 0.01 ~ "**",
                       p < 0.05 ~ "*"
                   ))
        
        wilcox_result <- data %>%
            semi_join(valid_groups[[analysis_type]], by = c("Tissue", "Type")) %>%
            group_by(Tissue) %>%
            filter(!is.na(Distance)) %>%
            filter(n_distinct(Type) >= 2) %>%
            wilcox_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(analysis = analysis_type,
                   test_type = "wilcox",
                   p.adj.signif = case_when(
                       p >= 0.05 ~ "ns",
                       p < 0.001 ~ "***",
                       p < 0.01 ~ "**",
                       p < 0.05 ~ "*"
                   ))
        
        return(list(t_test = t_test_result, wilcox = wilcox_result))
    }
    return(NULL)
}

#' Process Tissue-Specific Statistical Data
#'
#' @param file_path Path to RData file containing tissue statistics
#' @param stats_pattern Pattern to match statistics objects in loaded data
#' @param type_suffix Suffix to remove from data names to get type names
#'
#' @return List containing two data frames:
#'   \itemize{
#'     \item mean: Tissue-specific mean distances with columns Family, Type, Tissue, Distance
#'     \item median: Tissue-specific median distances with same columns
#'   }
#'
#' @importFrom dplyr select rename mutate bind_rows
#'
#' @examples
#' \dontrun{
#' stats <- process_tissue_statistics(
#'   "data/tissue_stats.RData",
#'   "^tissue_stats_",
#'   "_tissue_stats$"
#' )
#' }
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

#' @param df_median.dists_tissue Data frame of tissue-specific median distances
#' @param results_dir Directory path for saving output files
#' @param is_log2 Logical indicating if data is log2-transformed (default: FALSE)
#'
#' @return List containing:
#'   \itemize{
#'     \item t_test: List of mean and median t-test results
#'     \item wilcox: List of mean and median Wilcoxon test results
#'     \item summary: Combined test results data frame
#'   }
#'
#' @importFrom dplyr bind_rows mutate group_by summarise filter select
#' @importFrom purrr map
#' @importFrom utils write.csv
#'
#' @examples
#' \dontrun{
#' results <- perform_tissue_statistical_analysis(
#'   mean_data,
#'   median_data,
#'   "results/",
#'   is_log2 = FALSE
#' )
#' }
perform_tissue_statistical_analysis <- function(df_mean.dists_tissue, df_median.dists_tissue, results_dir, is_log2 = FALSE) {
    # Create output filename with appropriate suffix
    output_filename <- if(is_log2) {
        "statistical_tests_summary_tissue_log2.csv"
    } else {
        "statistical_tests_summary_tissue_regular.csv"
    }

    # Validate and prepare groups for testing
    valid_groups_tissue <- bind_rows(
        df_mean.dists_tissue %>% mutate(source = "mean"),
        df_median.dists_tissue %>% mutate(source = "median")
    ) %>%
        group_by(Tissue, Type, source) %>%
        summarise(n = n(), .groups = 'drop') %>%
        filter(n > 1) %>%
        split(.$source) %>%
        map(~select(.x, Tissue, Type))

    # Perform statistical tests
    test_results_tissue <- tryCatch({
        median_results <- perform_tissue_tests(df_median.dists_tissue, valid_groups_tissue, "median")
        mean_results <- perform_tissue_tests(df_mean.dists_tissue, valid_groups_tissue, "mean")
        
        if (!is.null(median_results) && !is.null(mean_results)) {
            test_summary_tissue <- bind_rows(
                median_results$t_test, mean_results$t_test,
                median_results$wilcox, mean_results$wilcox
            )
            write.csv(test_summary_tissue,
                     file.path(results_dir, output_filename),
                     row.names = FALSE)
            message(paste("Tissue-specific statistical tests summary exported to", output_filename))
        }
        
        list(
            t_test = list(
                median = if (!is.null(median_results)) median_results$t_test else NULL,
                mean = if (!is.null(mean_results)) mean_results$t_test else NULL
            ),
            wilcox = list(
                median = if (!is.null(median_results)) median_results$wilcox else NULL,
                mean = if (!is.null(mean_results)) mean_results$wilcox else NULL
            ),
            summary = if (exists("test_summary_tissue")) test_summary_tissue else NULL
        )
    }, error = function(e) {
        message("Error in tissue-specific statistical tests: ", e$message)
        return(NULL)
    })
    
    return(test_results_tissue)
}

# -------------------------------------------------------------------

#' Process Tissue-Specific Distance Data
#'
#' @param data_pattern Pattern to match distance data objects
#' @param loaded_objects List of loaded R objects
#'
#' @return Data frame containing processed tissue-specific distances with columns:
#'   \itemize{
#'     \item Family: Identifier for gene family
#'     \item Type: Type of distance measurement
#'     \item Tissue: Tissue type
#'     \item Distance: Distance value
#'   }
#'
#' @importFrom purrr map_df
#' @importFrom dplyr mutate select rename
#' @importFrom tidyr unnest_longer
#' @importFrom tibble enframe
#'
#' @examples
#' \dontrun{
#' tissue_distances <- process_tissue_distances(
#'   "^distance_data_",
#'   loaded_objects
#' )
#' }
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

#' Perform Statistical Analysis on Complete Tissue-Specific Dataset
#'
#' @param df_complete_dists Data frame containing complete tissue-specific distances
#' @param output_name Name of output CSV file
#' @param results_dir Directory path for saving output
#' @param is_log2 Logical indicating if data is log2-transformed (default: FALSE)
#'
#' @return List containing:
#'   \itemize{
#'     \item t_test: Results of tissue-specific t-tests
#'     \item wilcox: Results of tissue-specific Wilcoxon tests
#'   }
#' or NULL if analysis fails
#'
#' @importFrom dplyr group_by summarise filter select bind_rows
#' @importFrom utils write.csv
#'
#' @examples
#' \dontrun{
#' results <- perform_tissue_statistical_analysis_complete(
#'   tissue_data,
#'   "tissue_stats.csv",
#'   "results/"
#' )
#' }
perform_tissue_statistical_analysis_complete <- function(df_complete_dists, output_name, results_dir, is_log2 = FALSE) {
    valid_groups <- list(
        complete = df_complete_dists %>%
            group_by(Tissue, Type) %>%
            summarise(n = n(), .groups = 'drop') %>%
            filter(n > 1) %>%
            select(Tissue, Type)
    )

    test_results_tissue <- tryCatch({
        results <- perform_tissue_tests(df_complete_dists, valid_groups, "complete")
        
        if (!is.null(results)) {
            test_summary_tissue <- bind_rows(
                results$t_test,
                results$wilcox
            )
            write.csv(test_summary_tissue,
                     file.path(results_dir, output_name),
                     row.names = FALSE)
            message(paste("Tissue-specific statistical tests summary exported to", output_name))
        }
        
        results
    }, error = function(e) {
        message("Error in tissue-specific statistical tests: ", e$message)
        return(NULL)
    })

    return(test_results_tissue)
}

