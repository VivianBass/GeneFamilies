
library(testthat)
library(GeneFamilies)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

source("../../R/load_data_funks.R")

test_that("load_data_frame loads data with correct structure", {
    # Create a temporary file with test data
    temp_file <- tempfile()
    writeLines(
        "Family\tGene\tGene_species\tOrtholog\tOrtholog_species\nF1\tG1\tSp1\tO1\tSpO1\nF2\tG2\tSp2\tO2\tSpO2",
        temp_file
    )
    
    # Load the data using the function
    df <- load_data_frame(temp_file)
    
    # Check that the data frame has the correct column names
    expect_equal(colnames(df), c("Family", "Gene", "Gene_species", "Ortholog", "Ortholog_species"))
    
    # Check that the data frame has the correct number of rows
    expect_equal(nrow(df), 2)
    
    # Clean up the temporary file
    unlink(temp_file)
})


# test-load_data_frame.R
library(testthat)


test_that("load_data_frame loads valid file correctly", {
    temp_file <- tempfile()
    writeLines("Family\tGene\tGene_species\tOrtholog\tOrtholog_species\nF1\tG1\tSp1\tO1\tSpO1", temp_file)
    df <- load_data_frame(temp_file)
    expect_equal(colnames(df), c("Family", "Gene", "Gene_species", "Ortholog", "Ortholog_species"))
    expect_equal(nrow(df), 1)
    expect_equal(df$Family[1], "F1")
    unlink(temp_file)
})

test_that("load_data_frame handles missing columns", {
    temp_file <- tempfile()
    writeLines("Family\tGene\tGene_species\tOrtholog\nF1\tG1\tSp1\tO1", temp_file)
    expect_error(load_data_frame(temp_file), "more columns than column names") # Expected error due to missing column
    unlink(temp_file)
})

test_that("load_data_frame handles incorrect column names", {
    temp_file <- tempfile()
    writeLines("Family\tGene\tGene_species\tWrongName\tOrtholog_species\nF1\tG1\tSp1\tO1\tSpO1", temp_file)
    expect_error(load_data_frame(temp_file), "incorrect column names") # Expected error due to mismatched names
    unlink(temp_file)
})

test_that("load_data_frame loads empty file with headers", {
    temp_file <- tempfile()
    writeLines("Family\tGene\tGene_species\tOrtholog\tOrtholog_species\n", temp_file)
    df <- load_data_frame(temp_file)
    expect_equal(nrow(df), 0)
    expect_equal(colnames(df), c("Family", "Gene", "Gene_species", "Ortholog", "Ortholog_species"))
    unlink(temp_file)
})

test_that("load_data_frame handles non-tab delimited file", {
    temp_file <- tempfile()
    writeLines("Family,Gene,Gene_species,Ortholog,Ortholog_species\nF1,G1,Sp1,O1,SpO1", temp_file)
    df <- load_data_frame(temp_file)
    expect_false(all(colnames(df) == c("Family", "Gene", "Gene_species", "Ortholog", "Ortholog_species")))
    unlink(temp_file)
})

test_that("load_data_frame handles extra columns gracefully", {
    temp_file <- tempfile()
    writeLines("Family\tGene\tGene_species\tOrtholog\tOrtholog_species\tExtraColumn\nF1\tG1\tSp1\tO1\tSpO1\textra", temp_file)
    df <- load_data_frame(temp_file)
    expect_equal(ncol(df), 5) # Only 5 columns should be loaded
    unlink(temp_file)
})

test_that("load_data_frame handles file not found error", {
    expect_error(load_data_frame("non_existent_file.txt"), "cannot open file")
})
