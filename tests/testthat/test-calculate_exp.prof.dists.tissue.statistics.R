
library(testthat)
library(GeneFamilies)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

source("../../R/compute_funks.R")

test_that("calculate_exp.prof.dists.tissue.statistics computes mean and median per tissue for each family", {
    # Sample data with valid distance lists for tissues
    data <- list(
        Cluster1 = list(tissue1 = c(1, 2, 3), tissue2 = c(4, 5, 6)),
        Cluster2 = list(tissue1 = c(7, 8, 9), tissue2 = c(10, 11, 12))
    )
    
    # Apply the function
    result <- calculate_exp.prof.dists.tissue.statistics(data)
    
    # Expected means and medians for each cluster and tissue
    expected_means <- c(mean(c(1, 2, 3)), mean(c(4, 5, 6)), mean(c(7, 8, 9)), mean(c(10, 11, 12)))
    expected_medians <- c(median(c(1, 2, 3)), median(c(4, 5, 6)), median(c(7, 8, 9)), median(c(10, 11, 12)))
    
    # Check that result contains correct rows and values for means and medians
    expect_equal(result$Mean, expected_means)
    expect_equal(result$Median, expected_medians)
    expect_equal(result$Family, c("Cluster1", "Cluster1", "Cluster2", "Cluster2"))
    expect_equal(result$Tissue, c("tissue1", "tissue2", "tissue1", "tissue2"))
})

test_that("calculate_exp.prof.dists.tissue.statistics handles tissues with NA values", {
    # Sample data with NA values in tissue distance lists
    data <- list(
        Cluster1 = list(tissue1 = c(1, NA, 3), tissue2 = c(4, 5, NA)),
        Cluster2 = list(tissue1 = c(7, 8, 9), tissue2 = c(NA, 11, 12))
    )
    
    # Apply the function
    result <- calculate_exp.prof.dists.tissue.statistics(data)
    
    # Expected values excluding NA
    expected_means <- c(mean(c(1, 3), na.rm = TRUE), mean(c(4, 5), na.rm = TRUE), 
                        mean(c(7, 8, 9), na.rm = TRUE), mean(c(11, 12), na.rm = TRUE))
    expected_medians <- c(median(c(1, 3), na.rm = TRUE), median(c(4, 5), na.rm = TRUE), 
                          median(c(7, 8, 9), na.rm = TRUE), median(c(11, 12), na.rm = TRUE))
    
    # Check that result correctly computes mean and median excluding NAs
    expect_equal(result$Mean, expected_means)
    expect_equal(result$Median, expected_medians)
})

test_that("calculate_exp.prof.dists.tissue.statistics returns empty tibble for empty data", {
    # Empty data input
    data <- list()
    
    # Apply the function
    result <- calculate_exp.prof.dists.tissue.statistics(data)
    
    # Expect an empty tibble with columns Family, Tissue, Mean, and Median
    expect_equal(nrow(result), 0)
    expect_named(result, c("Family", "Tissue", "Mean", "Median"))
})

test_that("calculate_exp.prof.dists.tissue.statistics filters out infinite values", {
    # Sample data with Inf values in tissue lists
    data <- list(
        Cluster1 = list(tissue1 = c(1, Inf, 3), tissue2 = c(4, Inf, 6)),
        Cluster2 = list(tissue1 = c(7, 8, Inf), tissue2 = c(10, 11, Inf))
    )
    
    # Apply the function
    result <- calculate_exp.prof.dists.tissue.statistics(data)
    
    # Expected means and medians excluding Inf values
    expected_means <- c(mean(c(1, 3), na.rm = TRUE), mean(c(4, 6), na.rm = TRUE), 
                        mean(c(7, 8), na.rm = TRUE), mean(c(10, 11), na.rm = TRUE))
    expected_medians <- c(median(c(1, 3), na.rm = TRUE), median(c(4, 6), na.rm = TRUE), 
                          median(c(7, 8), na.rm = TRUE), median(c(10, 11), na.rm = TRUE))
    
    # Check that result correctly filters Inf values and computes mean and median
    expect_equal(result$Mean, expected_means)
    expect_equal(result$Median, expected_medians)
})

test_that("calculate_exp.prof.dists.tissue.statistics handles clusters with single-value tissues", {
    # Sample data with single-value distances per tissue
    data <- list(
        Cluster1 = list(tissue1 = c(5), tissue2 = c(10)),
        Cluster2 = list(tissue1 = c(15), tissue2 = c(20))
    )
    
    # Apply the function
    result <- calculate_exp.prof.dists.tissue.statistics(data)
    
    # Expected values (mean and median are the same for single-value lists)
    expected_means <- c(5, 10, 15, 20)
    expected_medians <- c(5, 10, 15, 20)
    
    # Check that result correctly handles single-value lists for mean and median
    expect_equal(result$Mean, expected_means)
    expect_equal(result$Median, expected_medians)
})
