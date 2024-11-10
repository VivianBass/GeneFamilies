
library(testthat)
library(GeneFamilies)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

source("../../R/compute_funks.R")

test_that("calculate_exp.prof.dists.statistics computes mean and median for each family", {
    # Sample data with valid distance matrices
    data <- list(
        Family1 = matrix(c(1, 2, 3, 4), nrow = 2),
        Family2 = matrix(c(5, 6, 7, 8), nrow = 2)
    )
    
    # Apply the function
    result <- calculate_exp.prof.dists.statistics(data)
    
    # Expected results
    expected_means <- c(mean(c(1, 2, 3, 4)), mean(c(5, 6, 7, 8)))
    expected_medians <- c(median(c(1, 2, 3, 4)), median(c(5, 6, 7, 8)))
    
    # Check that the result contains correct mean and median for each family
    expect_equal(result$Mean, expected_means)
    expect_equal(result$Median, expected_medians)
    expect_equal(result$Family, c("Family1", "Family2"))
})

test_that("calculate_exp.prof.dists.statistics handles matrices with NA values", {
    # Sample data with NA values in matrices
    data <- list(
        Family1 = matrix(c(1, NA, 3, 4), nrow = 2),
        Family2 = matrix(c(NA, 6, NA, 8), nrow = 2)
    )
    
    # Apply the function
    result <- calculate_exp.prof.dists.statistics(data)
    
    # Expected results with NA removed
    expected_means <- c(mean(c(1, 3, 4), na.rm = TRUE), mean(c(6, 8), na.rm = TRUE))
    expected_medians <- c(median(c(1, 3, 4), na.rm = TRUE), median(c(6, 8), na.rm = TRUE))
    
    # Check that the result correctly computes mean and median ignoring NAs
    expect_equal(result$Mean, expected_means)
    expect_equal(result$Median, expected_medians)
})

test_that("calculate_exp.prof.dists.statistics returns empty tibble for empty data", {
    # Empty data input
    data <- list()
    
    # Apply the function
    result <- calculate_exp.prof.dists.statistics(data)
    
    # Expect an empty tibble with columns Family, Mean, and Median
    expect_equal(nrow(result), 0)
    expect_named(result, c("Family", "Mean", "Median"))
})

test_that("calculate_exp.prof.dists.statistics filters out infinite values", {
    # Sample data with Inf values
    data <- list(
        Family1 = matrix(c(1, Inf, 3, 4), nrow = 2),
        Family2 = matrix(c(5, 6, Inf, 8), nrow = 2)
    )
    
    # Apply the function
    result <- calculate_exp.prof.dists.statistics(data)
    
    # Expected results with Inf removed
    expected_means <- c(mean(c(1, 3, 4), na.rm = TRUE), mean(c(5, 6, 8), na.rm = TRUE))
    expected_medians <- c(median(c(1, 3, 4), na.rm = TRUE), median(c(5, 6, 8), na.rm = TRUE))
    
    # Check that the result correctly filters out infinite values and computes mean and median
    expect_equal(result$Mean, expected_means)
    expect_equal(result$Median, expected_medians)
})

test_that("calculate_exp.prof.dists.statistics handles single-value matrices", {
    # Data with single values in matrices
    data <- list(
        Family1 = matrix(5, nrow = 1),
        Family2 = matrix(10, nrow = 1)
    )
    
    # Apply the function
    result <- calculate_exp.prof.dists.statistics(data)
    
    # Expected means and medians (both should be the single value in each case)
    expected_means <- c(5, 10)
    expected_medians <- c(5, 10)
    
    # Check that the result correctly computes mean and median for single-value matrices
    expect_equal(result$Mean, expected_means)
    expect_equal(result$Median, expected_medians)
})
