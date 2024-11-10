
library(testthat)
library(GeneFamilies)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

source("../../R/compute_funks.R")

test_that("exp.prof.dists calculates Euclidean distances for multiple genes", {
    # Define a sample expression profile with multiple genes across tissues
    expression.profiles <- data.frame(
        FBpp_ID = c("gene1", "gene2", "gene3"),
        tissue1 = c(1.2, 2.3, 3.1),
        tissue2 = c(2.5, 3.5, 4.1)
    )
    gene.accessions <- list("gene1", "gene2", "gene3")
    
    # Apply the function
    distances <- exp.prof.dists(gene.accessions, expression.profiles)
    
    # Check that distances vector is correct length for three genes (3 pairs)
    expect_length(distances, 3)
    expect_true(all(distances >= 0))  # Distances should be non-negative
})

test_that("exp.prof.dists returns NA if only one gene is provided", {
    # Define expression profiles with only one gene available
    expression.profiles <- data.frame(
        FBpp_ID = c("gene1"),
        tissue1 = c(1.2),
        tissue2 = c(2.5)
    )
    gene.accessions <- list("gene1")
    
    # Apply the function
    distances <- exp.prof.dists(gene.accessions, expression.profiles)
    
    # Expect NA as there's only one gene
    expect_true(is.na(distances))
})

test_that("exp.prof.dists handles missing columns gracefully", {
    # Define expression profiles missing one of the specified tissue columns
    expression.profiles <- data.frame(
        FBpp_ID = c("gene1", "gene2"),
        tissue1 = c(1.2, 2.3)
    )
    gene.accessions <- list("gene1", "gene2")
    
    # Apply the function with a non-existent column in tissues
    distances <- exp.prof.dists(gene.accessions, expression.profiles, tissues = c("tissue1", "tissue2"))
    
    # Expect an error due to missing column
    expect_error(distances, "undefined columns selected")
})

test_that("exp.prof.dists computes distances using specified distance method", {
    # Define a sample expression profile with genes and multiple tissues
    expression.profiles <- data.frame(
        FBpp_ID = c("gene1", "gene2", "gene3"),
        tissue1 = c(1.2, 2.3, 3.1),
        tissue2 = c(2.5, 3.5, 4.1)
    )
    gene.accessions <- list("gene1", "gene2", "gene3")
    
    # Apply the function with "manhattan" distance method
    distances <- exp.prof.dists(gene.accessions, expression.profiles, dist.method = "manhattan")
    
    # Check that distances vector is correct length for three genes
    expect_length(distances, 3)
})

test_that("exp.prof.dists returns empty vector if no genes match", {
    # Define an expression profile with no matching genes
    expression.profiles <- data.frame(
        FBpp_ID = c("gene4", "gene5"),
        tissue1 = c(1.0, 1.5),
        tissue2 = c(2.0, 2.5)
    )
    gene.accessions <- list("gene1", "gene2", "gene3")
    
    # Apply the function
    distances <- exp.prof.dists(gene.accessions, expression.profiles)
    
    # Expect NA as there are no matching genes
    expect_true(is.na(distances))
})
