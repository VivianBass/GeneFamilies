
library(testthat)
library(GeneFamilies)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

source("../../R/load_data_funks.R")

test_that("filter_v1 returns only matching genes", {
    # Define a sample df and expression_data with some overlapping genes
    df <- data.frame(
        Gene = c("gene1", "gene2", "gene3", "gene4"),
        Value = c(10, 20, 30, 40)
    )
    expression_data <- data.frame(
        FBpp_ID = c("gene2", "gene3", "gene5")
    )
    
    # Apply the function
    filtered_df <- filter_v1(df, expression_data)
    
    # Expect the result to only contain rows for "gene2" and "gene3"
    expect_equal(filtered_df$Gene, c("gene2", "gene3"))
    expect_equal(nrow(filtered_df), 2)
})

test_that("filter_v1 returns empty data frame if no genes match", {
    # Define a df and expression_data with no overlapping genes
    df <- data.frame(
        Gene = c("gene1", "gene2"),
        Value = c(10, 20)
    )
    expression_data <- data.frame(
        FBpp_ID = c("gene3", "gene4")
    )
    
    # Apply the function
    filtered_df <- filter_v1(df, expression_data)
    
    # Expect an empty data frame in this case
    expect_equal(nrow(filtered_df), 0)
})

test_that("filter_v1 returns all rows if all genes match", {
    # Define a df and expression_data with complete overlap
    df <- data.frame(
        Gene = c("gene1", "gene2"),
        Value = c(10, 20)
    )
    expression_data <- data.frame(
        FBpp_ID = c("gene1", "gene2")
    )
    
    # Apply the function
    filtered_df <- filter_v1(df, expression_data)
    
    # Expect the result to contain all rows
    expect_equal(nrow(filtered_df), nrow(df))
    expect_equal(filtered_df$Gene, df$Gene)
})

test_that("filter_v1 handles empty input data frame", {
    # Define an empty df and expression_data
    df <- data.frame(Gene = character(), Value = numeric())
    expression_data <- data.frame(FBpp_ID = c("gene1", "gene2"))
    
    # Apply the function
    filtered_df <- filter_v1(df, expression_data)
    
    # Expect an empty data frame in return
    expect_equal(nrow(filtered_df), 0)
})

test_that("filter_v1 handles empty expression_data frame", {
    # Define df with genes and an empty expression_data
    df <- data.frame(
        Gene = c("gene1", "gene2"),
        Value = c(10, 20)
    )
    expression_data <- data.frame(FBpp_ID = character())
    
    # Apply the function
    filtered_df <- filter_v1(df, expression_data)
    
    # Expect an empty data frame in return
    expect_equal(nrow(filtered_df), 0)
})
