
library(testthat)
library(GeneFamilies)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

source("../../R/load_data_funks.R")

# test-create_nested_list.R
test_that("create_nested_list creates nested list for Orthologs", {
    # Create a sample data frame
    df <- data.frame(
        Family = c("F1", "F1"),
        Gene = c("G1", "G2"),
        Gene_species = c("Sp1", "Sp1"),
        Ortholog = c("O1", "O2"),
        Ortholog_species = c("SpO1", "SpO2"),
        stringsAsFactors = FALSE
    )
    
    # Call the function with "Ortholog" header_type
    nested_list <- create_nested_list(df, "Ortholog")
    
    # Define the expected nested list structure
    expected <- list(
        F1 = list(
            "(Sp1, G1)" = list(SpO1 = "O1"),
            "(Sp1, G2)" = list(SpO2 = "O2")
        )
    )
    
    # Test if the output matches the expected structure
    expect_equal(nested_list, expected)
})

test_that("create_nested_list creates nested list for Paralogs", {
    # Create a sample data frame
    df <- data.frame(
        Family = c("F1", "F1"),
        Gene = c("G1", "G2"),
        Gene_species = c("Sp1", "Sp1"),
        Paralog = c("P1", "P2"),
        Paralog_species = c("SpP1", "SpP2"),
        stringsAsFactors = FALSE
    )
    
    # Call the function with "Paralog" header_type
    nested_list <- create_nested_list(df, "Paralog")
    
    # Define the expected nested list structure
    expected <- list(
        F1 = list(
            "(Sp1, G1)" = list(SpP1 = "P1"),
            "(Sp1, G2)" = list(SpP2 = "P2")
        )
    )
    
    # Test if the output matches the expected structure
    expect_equal(nested_list, expected)
})

test_that("create_nested_list handles an empty data frame", {
    # Create an empty data frame
    df <- data.frame(
        Family = character(),
        Gene = character(),
        Gene_species = character(),
        Ortholog = character(),
        Ortholog_species = character(),
        stringsAsFactors = FALSE
    )
    
    # Call the function with "Ortholog" header_type
    nested_list <- create_nested_list(df, "Ortholog")
    
    # Expect an empty list as the result
    expect_equal(nested_list, list())
})

test_that("create_nested_list creates separate nested lists for multiple families", {
    # Create a data frame with multiple families
    df <- data.frame(
        Family = c("F1", "F2", "F1", "F2"),
        Gene = c("G1", "G2", "G3", "G4"),
        Gene_species = c("Sp1", "Sp2", "Sp1", "Sp2"),
        Ortholog = c("O1", "O2", "O3", "O4"),
        Ortholog_species = c("SpO1", "SpO2", "SpO1", "SpO2"),
        stringsAsFactors = FALSE
    )
    
    # Call the function with "Ortholog" header_type
    nested_list <- create_nested_list(df, "Ortholog")
    
    # Define the expected nested list structure for multiple families
    expected <- list(
        F1 = list(
            "(Sp1, G1)" = list(SpO1 = "O1"),
            "(Sp1, G3)" = list(SpO1 = "O3")
        ),
        F2 = list(
            "(Sp2, G2)" = list(SpO2 = "O2"),
            "(Sp2, G4)" = list(SpO2 = "O4")
        )
    )
    
    # Test if the output matches the expected structure
    expect_equal(nested_list, expected)
})

test_that("create_nested_list throws error for missing required columns", {
    # Create a data frame with missing columns
    df <- data.frame(
        Family = c("F1", "F1"),
        Gene_species = c("Sp1", "Sp1"),
        Ortholog_species = c("SpO1", "SpO2"),
        stringsAsFactors = FALSE
    )
    
    # Expect an error due to missing columns
    expect_error(create_nested_list(df, "Ortholog"), "Required columns are missing.")
})

test_that("create_nested_list handles different header types", {
    # Create a sample data frame
    df <- data.frame(
        Family = c("F1", "F1"),
        Gene = c("G1", "G2"),
        Gene_species = c("Sp1", "Sp1"),
        Paralog = c("P1", "P2"),
        Paralog_species = c("SpP1", "SpP2"),
        stringsAsFactors = FALSE
    )
    
    # Call the function with "Paralog" header_type
    nested_list <- create_nested_list(df, "Paralog")
    
    # Define the expected structure for "Paralog" header type
    expected <- list(
        F1 = list(
            "(Sp1, G1)" = list(SpP1 = "P1"),
            "(Sp1, G2)" = list(SpP2 = "P2")
        )
    )
    
    # Test if the output matches the expected structure
    expect_equal(nested_list, expected)
})

test_that("create_nested_list handles duplicate genes within the same family", {
    # Create a sample data frame with duplicate genes
    df <- data.frame(
        Family = c("F1", "F1"),
        Gene = c("G1", "G1"),
        Gene_species = c("Sp1", "Sp1"),
        Ortholog = c("O1", "O2"),
        Ortholog_species = c("SpO1", "SpO2"),
        stringsAsFactors = FALSE
    )
    
    # Call the function with "Ortholog" header_type
    nested_list <- create_nested_list(df, "Ortholog")
    
    # Define expected structure, assuming duplicates are merged
    expected <- list(
        F1 = list(
            "(Sp1, G1)" = list(SpO1 = "O1", SpO2 = "O2")
        )
    )
    
    # Test if the output matches the expected structure
    expect_equal(nested_list, expected)
})

test_that("create_nested_list handles non-character data types", {
    # Create a sample data frame with non-character data types
    df <- data.frame(
        Family = c("F1", "F1"),
        Gene = as.factor(c("G1", "G2")),
        Gene_species = c("Sp1", "Sp1"),
        Ortholog = as.integer(c(1, 2)),
        Ortholog_species = c("SpO1", "SpO2"),
        stringsAsFactors = FALSE
    )
    
    # Call the function with "Ortholog" header_type
    nested_list <- create_nested_list(df, "Ortholog")
    
    # Define the expected structure, with integers converted to character
    expected <- list(
        F1 = list(
            "(Sp1, G1)" = list(SpO1 = "1"),
            "(Sp1, G2)" = list(SpO2 = "2")
        )
    )
    
    # Test if the output matches the expected structure
    expect_equal(nested_list, expected)
})

test_that("create_nested_list throws error for invalid header type", {
    # Create a sample data frame
    df <- data.frame(
        Family = c("F1"),
        Gene = c("G1"),
        Gene_species = c("Sp1"),
        Ortholog = c("O1"),
        Ortholog_species = c("SpO1"),
        stringsAsFactors = FALSE
    )
    
    # Test that an invalid header_type throws an error
    expect_error(create_nested_list(df, "InvalidType"), "Invalid header type provided.")
})
