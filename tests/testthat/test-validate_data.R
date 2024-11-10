
library(testthat)
library(GeneFamilies)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

source("../../R/compute_funks.R")

# Define test objects in the global environment for testing purposes
test_list_valid <- list(1, 2, 3)
test_vector_valid <- c(1, 2, 3)
test_list_na <- list(NA, NA, NA)
test_vector_na <- c(NA, NA, NA)
test_empty_list <- list()
test_invalid_type <- data.frame(a = 1:3)

# List of loaded objects' names as character vector for test purposes
loaded_objects <- c("test_list_valid", "test_vector_valid", "test_list_na", 
                    "test_vector_na", "test_empty_list", "test_invalid_type")

# Source the function if necessary
source("path/to/your/validate_data_function.R")

# Define pattern to match test object names ending in "_valid" or "_na"
pattern <- "_valid$|_na$"

test_that("validate_data correctly identifies valid lists and vectors", {
    # Expected output should only include test_list_valid and test_vector_valid
    valid_data_names <- validate_data(loaded_objects, pattern)
    expect_true("test_list_valid" %in% valid_data_names)
    expect_true("test_vector_valid" %in% valid_data_names)
})

test_that("validate_data excludes objects with only NA values", {
    # Expected output should exclude objects with only NA values
    valid_data_names <- validate_data(loaded_objects, pattern)
    expect_false("test_list_na" %in% valid_data_names)
    expect_false("test_vector_na" %in% valid_data_names)
})

test_that("validate_data excludes empty lists", {
    # Expected output should exclude empty lists
    valid_data_names <- validate_data(loaded_objects, pattern)
    expect_false("test_empty_list" %in% valid_data_names)
})

test_that("validate_data excludes invalid data types", {
    # Expected output should exclude data frames or non-list/vector types
    valid_data_names <- validate_data(loaded_objects, pattern)
    expect_false("test_invalid_type" %in% valid_data_names)
})

test_that("validate_data returns empty vector if no objects match pattern", {
    # Pattern that doesn't match any objects
    pattern_non_match <- "_nonexistent$"
    valid_data_names <- validate_data(loaded_objects, pattern_non_match)
    expect_equal(valid_data_names, character(0))
})

test_that("validate_data returns empty vector if no valid objects", {
    # Pattern that matches only invalid/empty objects
    pattern_invalid <- "_na$|_empty_list$"
    valid_data_names <- validate_data(loaded_objects, pattern_invalid)
    expect_equal(valid_data_names, character(0))
})
