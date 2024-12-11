library(testthat)
library(GeneFamilies)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

source("../../R/compute_funks.R")

test_that("calculate_exp.prof.dists.tissue.statistics computes mean and median per tissue for each family", {
    data <- list(
        Cluster1 = list(tissue1 = c(1, 2, 3), tissue2 = c(4, 5, 6)),
        Cluster2 = list(tissue1 = c(7, 8, 9), tissue2 = c(10, 11, 12))
    )
    
    result <- calculate_exp.prof.dists.tissue.statistics(data)
    
    expected_means <- setNames(c(mean(c(1, 2, 3)), mean(c(4, 5, 6)), 
                                mean(c(7, 8, 9)), mean(c(10, 11, 12))),
                              c("tissue1", "tissue2", "tissue1", "tissue2"))
    expected_medians <- setNames(c(median(c(1, 2, 3)), median(c(4, 5, 6)), 
                                  median(c(7, 8, 9)), median(c(10, 11, 12))),
                                c("tissue1", "tissue2", "tissue1", "tissue2"))
    
    expect_equal(result$Mean, expected_means)
    expect_equal(result$Median, expected_medians)
    expect_equal(result$Family, c("Cluster1", "Cluster1", "Cluster2", "Cluster2"))
    expect_equal(result$Tissue, c("tissue1", "tissue2", "tissue1", "tissue2"))
})

test_that("calculate_exp.prof.dists.tissue.statistics handles tissues with NA values", {
    data <- list(
        Cluster1 = list(tissue1 = c(1, NA, 3), tissue2 = c(4, 5, NA)),
        Cluster2 = list(tissue1 = c(7, 8, 9), tissue2 = c(NA, 11, 12))
    )
    
    result <- calculate_exp.prof.dists.tissue.statistics(data)
    
    expected_means <- setNames(c(mean(c(1, 3), na.rm = TRUE), mean(c(4, 5), na.rm = TRUE),
                                mean(c(7, 8, 9), na.rm = TRUE), mean(c(11, 12), na.rm = TRUE)),
                              c("tissue1", "tissue2", "tissue1", "tissue2"))
    expected_medians <- setNames(c(median(c(1, 3), na.rm = TRUE), median(c(4, 5), na.rm = TRUE),
                                  median(c(7, 8, 9), na.rm = TRUE), median(c(11, 12), na.rm = TRUE)),
                                c("tissue1", "tissue2", "tissue1", "tissue2"))
    
    expect_equal(result$Mean, expected_means)
    expect_equal(result$Median, expected_medians)
})

test_that("calculate_exp.prof.dists.tissue.statistics returns empty tibble for empty data", {
    data <- list()
    result <- calculate_exp.prof.dists.tissue.statistics(data)
    expect_equal(nrow(result), 0)
    expect_named(result, c("Family", "Tissue", "Mean", "Median"))
})

test_that("calculate_exp.prof.dists.tissue.statistics filters out infinite values", {
    data <- list(
        Cluster1 = list(tissue1 = c(1, Inf, 3), tissue2 = c(4, Inf, 6)),
        Cluster2 = list(tissue1 = c(7, 8, Inf), tissue2 = c(10, 11, Inf))
    )
    
    result <- calculate_exp.prof.dists.tissue.statistics(data)
    
    expected_means <- setNames(c(mean(c(1, 3), na.rm = TRUE), mean(c(4, 6), na.rm = TRUE),
                                mean(c(7, 8), na.rm = TRUE), mean(c(10, 11), na.rm = TRUE)),
                              c("tissue1", "tissue2", "tissue1", "tissue2"))
    expected_medians <- setNames(c(median(c(1, 3), na.rm = TRUE), median(c(4, 6), na.rm = TRUE),
                                  median(c(7, 8), na.rm = TRUE), median(c(10, 11), na.rm = TRUE)),
                                c("tissue1", "tissue2", "tissue1", "tissue2"))
    
    expect_equal(result$Mean, expected_means)
    expect_equal(result$Median, expected_medians)
})

test_that("calculate_exp.prof.dists.tissue.statistics handles clusters with single-value tissues", {
    data <- list(
        Cluster1 = list(tissue1 = c(5), tissue2 = c(10)),
        Cluster2 = list(tissue1 = c(15), tissue2 = c(20))
    )
    
    result <- calculate_exp.prof.dists.tissue.statistics(data)
    
    expected_means <- setNames(c(5, 10, 15, 20), c("tissue1", "tissue2", "tissue1", "tissue2"))
    expected_medians <- setNames(c(5, 10, 15, 20), c("tissue1", "tissue2", "tissue1", "tissue2"))
    
    expect_equal(result$Mean, expected_means)
    expect_equal(result$Median, expected_medians)
})
