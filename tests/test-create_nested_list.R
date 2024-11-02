

library(testthat)
library(GeneFamilies)

test_check("GeneFamilies")

# test-create_nested_list.R
library(testthat)
library(dplyr)

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
