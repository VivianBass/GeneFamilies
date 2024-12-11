library(testthat)
library(GeneFamilies)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

source("../../R/load_data_funks.R")

# Test for Ortholog data
test_that("load_data_frame loads ortholog data correctly", {
    temp_file <- tempfile()
    writeLines(
        "Family\tGene\tGene_species\tOrtholog\tOrtholog_species\nF1\tG1\tSp1\tO1\tSpO1\nF2\tG2\tSp2\tO2\tSpO2",
        temp_file
    )
    df <- load_data_frame(temp_file)
    expect_equal(colnames(df), c("Family", "Gene", "Gene_species", "Ortholog", "Ortholog_species"))
    expect_equal(nrow(df), 2)
    unlink(temp_file)
})

# Test for Paralog data
test_that("load_data_frame loads paralog data correctly", {
    temp_file <- tempfile()
    writeLines(
        "Family\tGene\tGene_species\tParalog\tParalog_species\nF1\tG1\tSp1\tP1\tSpP1\nF2\tG2\tSp2\tP2\tSpP2",
        temp_file
    )
    df <- load_data_frame(temp_file)
    expect_equal(colnames(df), c("Family", "Gene", "Gene_species", "Paralog", "Paralog_species"))
    expect_equal(nrow(df), 2)
    unlink(temp_file)
})

test_that("load_data_frame handles missing columns for both types", {
    # Test Ortholog missing columns
    temp_file1 <- tempfile()
    writeLines("Family\tGene\tGene_species\tOrtholog\nF1\tG1\tSp1\tO1", temp_file1)
    expect_error(load_data_frame(temp_file1), "incorrect column names")
    unlink(temp_file1)
    
    # Test Paralog missing columns
    temp_file2 <- tempfile()
    writeLines("Family\tGene\tGene_species\tParalog\nF1\tG1\tSp1\tP1", temp_file2)
    expect_error(load_data_frame(temp_file2), "incorrect column names")
    unlink(temp_file2)
})

test_that("load_data_frame handles incorrect column names for both types", {
    temp_file <- tempfile()
    writeLines("Family\tGene\tGene_species\tWrongName\tWrong_species\nF1\tG1\tSp1\tX1\tSpX1", temp_file)
    expect_error(load_data_frame(temp_file), "incorrect column names")
    unlink(temp_file)
})

test_that("load_data_frame loads empty files with valid headers", {
    # Test empty Ortholog file
    temp_file1 <- tempfile()
    writeLines("Family\tGene\tGene_species\tOrtholog\tOrtholog_species\n", temp_file1)
    df1 <- load_data_frame(temp_file1)
    expect_equal(nrow(df1), 0)
    expect_equal(colnames(df1), c("Family", "Gene", "Gene_species", "Ortholog", "Ortholog_species"))
    unlink(temp_file1)
    
    # Test empty Paralog file
    temp_file2 <- tempfile()
    writeLines("Family\tGene\tGene_species\tParalog\tParalog_species\n", temp_file2)
    df2 <- load_data_frame(temp_file2)
    expect_equal(nrow(df2), 0)
    expect_equal(colnames(df2), c("Family", "Gene", "Gene_species", "Paralog", "Paralog_species"))
    unlink(temp_file2)
})

test_that("load_data_frame handles extra columns gracefully for both types", {
    # Test Ortholog with extra column
    temp_file1 <- tempfile()
    writeLines("Family\tGene\tGene_species\tOrtholog\tOrtholog_species\tExtra\nF1\tG1\tSp1\tO1\tSpO1\textra", temp_file1)
    df1 <- load_data_frame(temp_file1)
    expect_equal(ncol(df1), 5)
    unlink(temp_file1)
    
    # Test Paralog with extra column
    temp_file2 <- tempfile()
    writeLines("Family\tGene\tGene_species\tParalog\tParalog_species\tExtra\nF1\tG1\tSp1\tP1\tSpP1\textra", temp_file2)
    df2 <- load_data_frame(temp_file2)
    expect_equal(ncol(df2), 5)
    unlink(temp_file2)
})
