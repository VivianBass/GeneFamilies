library(testthat)
library(GeneFamilies)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

source("../../R/angles_funks.R")

# Define test cases
test_that("Validates and retains non-empty data frames", {
  df_list <- list(
    df1 = data.frame(a = 1:3, b = 4:6),
    df2 = data.frame(x = c(10, 20), y = c(30, 40))
  )
  result <- validate_angle_dataframes(df_list)
  
  expect_equal(length(result), 2)
  expect_named(result, c("df1", "df2"))
  expect_s3_class(result$df1, "data.frame")
  expect_s3_class(result$df2, "data.frame")
})

test_that("Excludes null data frames", {
  df_list <- list(
    df1 = data.frame(a = 1:3, b = 4:6),
    df2 = NULL
  )
  expect_message(result <- validate_angle_dataframes(df_list), "Warning: df2 dataframe is empty or null")
  expect_equal(length(result), 1)
  expect_named(result, "df1")
})

test_that("Excludes empty data frames", {
  df_list <- list(
    df1 = data.frame(a = 1:3, b = 4:6),
    df2 = data.frame()
  )
  expect_message(result <- validate_angle_dataframes(df_list), "Warning: df2 dataframe is empty or null")
  expect_equal(length(result), 1)
  expect_named(result, "df1")
})

test_that("Handles a list with only null or empty data frames", {
  df_list <- list(
    df1 = NULL,
    df2 = data.frame()
  )
  
  messages <- testthat::capture_messages({
    result <- validate_angle_dataframes(df_list)
  })
  
  expect_true(any(grepl("Warning: df1 dataframe is empty or null", messages)))
  expect_true(any(grepl("Warning: df2 dataframe is empty or null", messages)))
  expect_equal(length(result), 0)
})

test_that("Processes a list with mixed valid and invalid entries", {
  df_list <- list(
    df1 = data.frame(a = 1:3, b = 4:6),
    df2 = NULL,
    df3 = data.frame(),
    df4 = data.frame(x = 1:2, y = 3:4)
  )
  
  messages <- testthat::capture_messages({
    result <- validate_angle_dataframes(df_list)
  })
  
  expect_true(any(grepl("Warning: df2 dataframe is empty or null", messages)))
  expect_true(any(grepl("Warning: df3 dataframe is empty or null", messages)))
  expect_equal(length(result), 2)
  expect_named(result, c("df1", "df4"))
})



test_that("Handles an empty input list", {
  df_list <- list()
  result <- validate_angle_dataframes(df_list)
  
  expect_equal(length(result), 0)
  expect_s3_class(result, "list")
  expect_equal(names(result), NULL)
})
