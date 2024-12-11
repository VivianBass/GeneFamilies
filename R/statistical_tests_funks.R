
#' Perform Statistical Tests on Distance Data
#' 
#' FUNCTION FOR ONE DIRECTIONAL TESTS
#'
#' @param data Data frame with columns 'Type' and 'Distance'
#' @param valid_groups List of valid groups for comparison
#' @param analysis_type String indicating "mean" or "median" analysis
#' @return List containing t-test and Wilcoxon test results, or NULL
#'
#' @importFrom dplyr filter mutate case_when
#' @importFrom rstatix t_test wilcox_test adjust_pvalue
#' Function is turned OFF by suffix _X
perform_tests_X <- function(data, valid_groups, analysis_type) {
    if (length(valid_groups[[analysis_type]]) >= 2) {

        t_test_result <- data %>%
            filter(Type %in% valid_groups[[analysis_type]]) %>%
            t_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(analysis = analysis_type, 
                   test_type = "t-test",
                   p.adj.signif = case_when(
                       p.adj >= 0.05 ~ "ns",
                       p.adj < 0.001 ~ "***",
                       p.adj < 0.01 ~ "**",
                       p.adj < 0.05 ~ "*"
                   ))
        
        wilcox_result <- data %>%
            filter(Type %in% valid_groups[[analysis_type]]) %>%
            wilcox_test(Distance ~ Type, alternative = "greater") %>%
            adjust_pvalue(method = "BH") %>%
            mutate(analysis = analysis_type, 
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


#' Perform Statistical Tests on Distance Data
#' 
#' FUNCTION FOR TWO DIRECTIONAL TESTS (group1 vs group2 and group2 vs group1)
#'
#' Performs t-tests and Wilcoxon tests on all possible group pairs in both directions
#' (group1 vs group2 and group2 vs group1).
#'
#' @param data Data frame with columns 'Type' and 'Distance'
#' @param valid_groups List of valid groups for analysis
#' @param analysis_type String specifying analysis type ("mean", "median", "complete")
#'
#' @return List containing t-test and Wilcoxon test results, or NULL if no valid results
#'
#' @details
#' - Requires ≥2 observations per group
#' - Adjusts p-values using BH method
#' - Orders results by test type
#'
#' @importFrom dplyr filter mutate arrange bind_rows
#' @importFrom rstatix t_test wilcox_test adjust_pvalue
#'
#' @examples
#' data <- data.frame(Type = c("A", "A", "B", "B"), Distance = c(1.2, 1.5, 2.3, 2.1))
#' results <- perform_tests(data, list(mean = c("A", "B")), "mean")
#'
#' @export
perform_tests <- function(data, valid_groups, analysis_type) {
    if (length(valid_groups[[analysis_type]]) >= 2) {
        # Create all possible combinations of groups in both directions
        groups <- valid_groups[[analysis_type]]
        group_pairs <- expand.grid(group1 = groups, group2 = groups)
        # Remove self-comparisons
        group_pairs <- group_pairs[group_pairs$group1 != group_pairs$group2,]
        
        # Initialize empty results
        t_test_results <- list()
        wilcox_results <- list()
        
        # Loop through each pair
        for(i in 1:nrow(group_pairs)) {
            g1 <- group_pairs$group1[i]
            g2 <- group_pairs$group2[i]
            
            # Get data for current pair
            test_data <- data %>% filter(Type %in% c(g1, g2))
            
            # Get counts for each group
            g1_count <- sum(test_data$Type == g1)
            g2_count <- sum(test_data$Type == g2)
            
            # Only perform tests if both groups have sufficient observations
            if(g1_count >= 2 && g2_count >= 2) {
                # T-test
                tryCatch({
                    t_result <- t_test(test_data, Distance ~ Type, alternative = "greater") %>%
                        adjust_pvalue(method = "BH") %>%
                        mutate(
                            analysis = analysis_type,
                            test_type = "t-test",
                            group1 = g1,
                            group2 = g2,
                            p.adj.signif = case_when(
                                p.adj >= 0.05 ~ "ns",
                                p.adj < 0.001 ~ "***",
                                p.adj < 0.01 ~ "**",
                                p.adj < 0.05 ~ "*"
                            )
                        )
                    t_test_results[[i]] <- t_result
                }, error = function(e) {
                    message(sprintf("T-test failed for %s vs %s: %s", g1, g2, e$message))
                })
                
                # Wilcoxon test
                tryCatch({
                    w_result <- wilcox_test(test_data, Distance ~ Type, alternative = "greater") %>%
                        adjust_pvalue(method = "BH") %>%
                        mutate(
                            analysis = analysis_type,
                            test_type = "wilcox",
                            group1 = g1,
                            group2 = g2,
                            p.adj.signif = case_when(
                                p.adj >= 0.05 ~ "ns",
                                p.adj < 0.001 ~ "***",
                                p.adj < 0.01 ~ "**",
                                p.adj < 0.05 ~ "*"
                            )
                        )
                    wilcox_results[[i]] <- w_result
                }, error = function(e) {
                    message(sprintf("Wilcoxon test failed for %s vs %s: %s", g1, g2, e$message))
                })
            } else {
                message(sprintf("Skipping comparison %s (n=%d) vs %s (n=%d): insufficient observations", 
                              g1, g1_count, g2, g2_count))
            }
        }
        
        # Combine results and order by test_type
        t_test_result <- bind_rows(t_test_results) %>% arrange(test_type)
        wilcox_result <- bind_rows(wilcox_results) %>% arrange(test_type)
        
        if(nrow(t_test_result) == 0 && nrow(wilcox_result) == 0) {
            message("No valid test results were produced")
            return(NULL)
        }
        
        return(list(t_test = t_test_result, wilcox = wilcox_result))
    }
    return(NULL)
}

# ---------------------------------------------------------------------------------

#' Process Mean and Median Distance Data
#'
#' @param input_file Path to input RData file
#' @param stats_pattern Pattern to match statistics objects
#' @param output_file Name for output file
#' @param output_data_dir Directory for output
#' @return List containing processed mean and median data frames
#'
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
process_distances <- function(input_file, stats_pattern, output_file, output_data_dir) {
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

#' Perform Statistical Tests on Mean and Median Data
#'
#' @param data_list List containing mean and median data frames
#' @param output_file Name for output CSV file
#' @return List of test results for mean and median data
#'
#' @importFrom dplyr group_by summarise filter pull bind_rows
#' @importFrom utils write.csv
perform_statistical_tests <- function(data_list, output_file) {
    # Extract mean and median dataframes from the list
    df_mean <- data_list$mean
    df_median <- data_list$median
    
    # Check for valid groups in mean data
    valid_groups <- df_mean %>%
        group_by(Type) %>%
        summarise(n = n(), .groups = 'drop') %>%
        filter(n > 1) %>%
        pull(Type)
    
    # Perform tests
    test_results <- tryCatch({
        # Perform tests for mean distances
        mean_tests <- perform_tests(df_mean, list(mean = valid_groups), "mean")
        
        # Perform tests for median distances
        median_tests <- perform_tests(df_median, list(median = valid_groups), "median")
        
        if (!is.null(mean_tests) && !is.null(median_tests)) {
            test_summary <- bind_rows(
                mean_tests$t_test,
                mean_tests$wilcox,
                median_tests$t_test,
                median_tests$wilcox
            )
            
            write.csv(test_summary,
                     file.path(results_dir, output_file),
                     row.names = FALSE)
            message("statistical tests summary exported to CSV")
        }
        
        list(
            mean_t_test = if(!is.null(mean_tests)) mean_tests$t_test else NULL,
            mean_wilcox = if(!is.null(mean_tests)) mean_tests$wilcox else NULL,
            median_t_test = if(!is.null(median_tests)) median_tests$t_test else NULL,
            median_wilcox = if(!is.null(median_tests)) median_tests$wilcox else NULL,
            summary = if(exists("test_summary")) test_summary else NULL
        )
    }, error = function(e) {
        message("Error in statistical tests: ", e$message)
        return(NULL)
    })
    
    return(test_results)
}

# ---------------------------------------------------------------------------------

#' Process Complete Distance Data
#'
#' @param data_pattern Pattern to match distance data objects
#' @param loaded_objects List of loaded R objects
#' @return Data frame of processed distances
#'
#' @importFrom purrr map_df
#' @importFrom tibble tibble
process_complete_distances <- function(data_pattern, loaded_objects) {
    data_names <- loaded_objects[grepl(data_pattern, loaded_objects)]
    
    data_names %>%
        map_df(function(data_name) {
            current_data <- get(data_name)
            type_name <- sub(data_pattern, "", data_name)
            
            # Check if we're handling angles or euclidean distances
            if (grepl("angles", data_name)) {
                tibble(
                    Type = type_name,
                    Angle = unlist(current_data)
                )
            } else {
                tibble(
                    Type = type_name,
                    Distance = unlist(current_data)
                )
            }
        })
}

#' Perform Statistical Tests on Complete Distance Data
#'
#' This function performs t-tests and Wilcoxon tests on a complete set of distance data, adjusting p-values and exporting results to a CSV file.
#'
#' @param df_complete A data frame containing the complete distance data with columns for Type and either Distance or Angle.
#' @param output_file A string specifying the name of the output CSV file where the test results will be saved.
#' @return A list containing the results of the t-tests and Wilcoxon tests, as well as a summary of the results.
#'
#' @importFrom dplyr group_by summarise filter pull rename
#' @importFrom rstatix t_test wilcox_test adjust_pvalue
#' @importFrom utils write.csv
#' @examples
#' # Example usage
#' df_complete <- data.frame(
#'   Type = c("A", "A", "B", "B"),
#'   Distance = c(1.2, 1.5, 2.3, 2.1)
#' )
#' results <- perform_statistical_tests_complete(df_complete, "output.csv")
#' @export
perform_statistical_tests_complete <- function(df_complete, output_file) {
    # Check for valid groups
    valid_groups <- df_complete %>%
        group_by(Type) %>%
        summarise(n = n(), .groups = 'drop') %>%
        filter(n > 1) %>%
        pull(Type)
    
    # Determine which metric to use based on column names
    metric_column <- if("Angle" %in% names(df_complete)) "Angle" else "Distance"
    
    # If metric is Angle, temporarily rename it to Distance for perform_tests
    if(metric_column == "Angle") {
        df_complete <- df_complete %>% 
            rename(Distance = Angle)
    }
    
    # Perform tests
    test_results <- tryCatch({
        # Create a list with valid groups for perform_tests
        valid_groups_list <- list(complete = valid_groups)
        
        # Use perform_tests
        test_results <- perform_tests(df_complete, valid_groups_list, "complete")
        
        if (!is.null(test_results)) {
            test_summary <- bind_rows(
                test_results$t_test,
                test_results$wilcox
            )
            
            write.csv(test_summary,
                     file.path(results_dir, output_file),
                     row.names = FALSE)
            message("Statistical tests summary exported to CSV")
        }
        
        list(
            t_test = if(!is.null(test_results)) test_results$t_test else NULL,
            wilcox = if(!is.null(test_results)) test_results$wilcox else NULL,
            summary = if(exists("test_summary")) test_summary else NULL
        )
    }, error = function(e) {
        message("Error in statistical tests: ", e$message)
        return(NULL)
    })
    
    return(test_results)
}