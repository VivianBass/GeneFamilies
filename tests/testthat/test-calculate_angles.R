
library(testthat)
library(GeneFamilies)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)


source("../../R/angles_funks.R")

# Load the function definitions
source("../../R/expression_funks.R")

# Test Cases for calculate_angles
test_that("calculate_angles computes angles correctly for valid input", {
  genes <- list(group1 = c("gene1", "gene2"))
  rna.seq.exp.profils <- data.frame(
    FBpp_ID = c("gene1", "gene2", "gene3"),
    tissue1 = c(1.0, 2.0, 3.0),
    tissue2 = c(1.0, 2.0, 3.0)
  )
  tissues <- c("tissue1", "tissue2")

  result <- calculate_angles(genes, rna.seq.exp.profils, tissues)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)  # Both genes should be included
  expect_named(result, c("FBpp_ID", "angle.diag"))
  expect_false(any(is.na(result$angle.diag)))
})

test_that("calculate_angles handles empty gene groups", {
  genes <- list()
  rna.seq.exp.profils <- data.frame(
    FBpp_ID = c("gene1", "gene2"),
    tissue1 = c(1.0, 2.0),
    tissue2 = c(1.0, 2.0)
  )
  tissues <- c("tissue1", "tissue2")
  
  result <- calculate_angles(genes, rna.seq.exp.profils, tissues)
  expect_equal(nrow(result), 0)
  expect_warning(calculate_angles(genes, rna.seq.exp.profils, tissues))
})

test_that("calculate_angles handles mismatched gene identifiers", {
  genes <- list(group1 = c("geneX", "geneY"))
  rna.seq.exp.profils <- data.frame(
    FBpp_ID = c("gene1", "gene2"),
    tissue1 = c(1.0, 2.0),
    tissue2 = c(1.0, 2.0)
  )
  tissues <- c("tissue1", "tissue2")
  
  expect_warning(result <- calculate_angles(genes, rna.seq.exp.profils, tissues))
  expect_equal(nrow(result), 0)
})
